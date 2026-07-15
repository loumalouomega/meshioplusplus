#pragma once
//
// NDArray: a minimal typed, n-dimensional, contiguous (row-major) array.
//
// It is the storage primitive used by meshioplusplus::Mesh for points, cell
// connectivity, and all data fields. An NDArray either *owns* its buffer
// (the common case, e.g. data produced by a reader) or is a non-owning
// *view* over external memory (used to wrap numpy buffers zero-copy when
// writing). The binding layer converts between NDArray and numpy.

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <memory>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

namespace meshioplusplus {

namespace detail {
// Allocator that leaves elements *default*-initialized rather than
// value-initialized. For a trivial type like std::byte that means the buffer
// is left uninitialized instead of zero-filled. NDArray uses it so a buffer it
// is about to fully overwrite (reader outputs, reconstruction blocks) can skip
// the zero-fill memset — which, for a fresh large allocation, is an entire
// extra cold pass over just-faulted pages (numpy's calloc-backed arrays skip
// it too). std::vector stays copyable, unlike a unique_ptr buffer.
template <class T>
struct no_init_allocator {
    using value_type = T;
    no_init_allocator() = default;
    template <class U>
    no_init_allocator(const no_init_allocator<U>&) noexcept {}
    template <class U>
    struct rebind {
        using other = no_init_allocator<U>;
    };
    T* allocate(std::size_t n) { return std::allocator<T>{}.allocate(n); }
    void deallocate(T* p, std::size_t n) { std::allocator<T>{}.deallocate(p, n); }
    // Default-init (no zeroing) for the no-arg case resize() uses; forward
    // everything else so the vector still behaves normally.
    template <class U>
    void construct(U* p) noexcept(std::is_nothrow_default_constructible_v<U>) {
        ::new (static_cast<void*>(p)) U;
    }
    template <class U, class... Args>
    void construct(U* p, Args&&... args) {
        ::new (static_cast<void*>(p)) U(std::forward<Args>(args)...);
    }
    template <class U>
    bool operator==(const no_init_allocator<U>&) const noexcept {
        return true;
    }
    template <class U>
    bool operator!=(const no_init_allocator<U>&) const noexcept {
        return false;
    }
};
}  // namespace detail

enum class DType {
    Float32,
    Float64,
    Int8,
    Int16,
    Int32,
    Int64,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
};

inline std::size_t dtype_size(DType dt) {
    switch (dt) {
        case DType::Float32: return 4;
        case DType::Float64: return 8;
        case DType::Int8:
        case DType::UInt8: return 1;
        case DType::Int16:
        case DType::UInt16: return 2;
        case DType::Int32:
        case DType::UInt32: return 4;
        case DType::Int64:
        case DType::UInt64: return 8;
    }
    return 0;
}

// numpy dtype string (kind + itemsize), e.g. "f8", "i4".
inline const char* dtype_numpy_str(DType dt) {
    switch (dt) {
        case DType::Float32: return "f4";
        case DType::Float64: return "f8";
        case DType::Int8: return "i1";
        case DType::Int16: return "i2";
        case DType::Int32: return "i4";
        case DType::Int64: return "i8";
        case DType::UInt8: return "u1";
        case DType::UInt16: return "u2";
        case DType::UInt32: return "u4";
        case DType::UInt64: return "u8";
    }
    return "f8";
}

class NDArray {
public:
    NDArray() = default;

    // Owning array, zero-initialised buffer of the right size.
    NDArray(DType dt, std::vector<std::size_t> shape)
        : dtype_(dt), shape_(std::move(shape)) {
        const std::size_t nb = nbytes();
        owned_.resize(nb);  // uninitialised (no_init_allocator)
        std::memset(owned_.data(), 0, nb);  // explicit zero-fill
    }

    // Owning array whose buffer is left *uninitialised* — only for callers that
    // immediately overwrite every byte (reader outputs, reconstruction blocks).
    static NDArray uninit(DType dt, std::vector<std::size_t> shape) {
        NDArray a;
        a.dtype_ = dt;
        a.shape_ = std::move(shape);
        a.owned_.resize(a.nbytes());  // no memset
        return a;
    }

    // Non-owning view over external row-major memory (caller keeps it alive).
    static NDArray make_view(DType dt, std::vector<std::size_t> shape, std::byte* ptr) {
        NDArray a;
        a.dtype_ = dt;
        a.shape_ = std::move(shape);
        a.view_ = ptr;
        return a;
    }

    DType dtype() const { return dtype_; }
    const std::vector<std::size_t>& shape() const { return shape_; }
    std::size_t ndim() const { return shape_.size(); }
    bool is_view() const { return view_ != nullptr; }

    std::size_t size() const {
        if (shape_.empty()) return 0;
        return std::accumulate(shape_.begin(), shape_.end(), std::size_t{1},
                               std::multiplies<std::size_t>());
    }
    std::size_t nbytes() const { return size() * dtype_size(dtype_); }

    std::byte* data() { return view_ ? view_ : owned_.data(); }
    const std::byte* data() const { return view_ ? view_ : owned_.data(); }

    // Change the logical shape without touching the buffer (sizes must match).
    void reshape(std::vector<std::size_t> new_shape) {
        std::size_t n = new_shape.empty()
                            ? 0
                            : std::accumulate(new_shape.begin(), new_shape.end(),
                                              std::size_t{1}, std::multiplies<std::size_t>());
        if (n != size()) return;  // ignore inconsistent reshape
        shape_ = std::move(new_shape);
    }

    // Turn a view into an owning copy (no-op if already owning). Used before
    // handing a buffer's lifetime to Python via a capsule.
    void make_owned() {
        if (view_ == nullptr) return;
        const std::size_t nb = nbytes();
        ByteBuf buf;
        buf.resize(nb);  // uninitialised; fully overwritten by the memcpy below
        std::memcpy(buf.data(), view_, nb);
        owned_ = std::move(buf);
        view_ = nullptr;
    }

    template <typename T>
    T* as() { return reinterpret_cast<T*>(data()); }
    template <typename T>
    const T* as() const { return reinterpret_cast<const T*>(data()); }

private:
    using ByteBuf = std::vector<std::byte, detail::no_init_allocator<std::byte>>;
    DType dtype_ = DType::Float64;
    std::vector<std::size_t> shape_;
    ByteBuf owned_;
    std::byte* view_ = nullptr;
};

}  // namespace meshioplusplus
