//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████                           
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███    
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███    
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░ 
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███    
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░     
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░                            
//                                                                                                  
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file ndarray.hpp
 * @brief `NDArray`: a minimal typed, n-dimensional, contiguous (row-major)
 * array — the storage primitive of `meshioplusplus::Mesh`.
 *
 * `NDArray` is used for points, cell connectivity, and every point/cell/field
 * data array. It either *owns* its buffer (the common case: data produced by
 * a reader) or is a non-owning *view* over externally-owned memory (used to
 * wrap a numpy buffer zero-copy on the write path — see `py_to_mesh` in
 * `bindings/np_conversions.hpp`). The binding layer converts between
 * `NDArray` and numpy at the I/O boundary: owning buffers are moved into a
 * capsule backing a writeable numpy array on read, and numpy buffers are
 * wrapped as views (no copy) on write. `dtype()` records the element type
 * with an internal `DType` enum rather than a template parameter, so
 * `NDArray` can be stored uniformly (e.g. in `Mesh::cell_data`) regardless of
 * the numpy dtype it came from; `as<T>()` reinterprets the raw buffer as `T`
 * once the caller has determined (or asserted) the appropriate type.
 */

// System includes
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
/**
 * @brief Allocator that leaves elements *default*-initialized rather than
 * value-initialized.
 *
 * For a trivial type like `std::byte` that means the buffer is left
 * uninitialized instead of zero-filled. `NDArray` uses this (via `ByteBuf`)
 * so a buffer it is about to fully overwrite (reader outputs, reconstruction
 * blocks — see `NDArray::uninit`) can skip the zero-fill `memset`, which for
 * a fresh large allocation is an entire extra cold pass over just-faulted
 * pages (numpy's `calloc`-backed arrays skip it too, for the same reason).
 * `std::vector` with this allocator stays copyable/movable like a normal
 * vector, unlike a raw `unique_ptr` buffer, so `NDArray` can keep value
 * semantics.
 *
 * @tparam T The element type being allocated (used as `std::byte` here).
 */
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

/**
 * @brief Scalar element type of an `NDArray`, mirroring the numpy dtypes the
 * binding layer converts to/from.
 */
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

/**
 * @brief Size in bytes of one element of the given dtype.
 * @param dt The dtype to query.
 * @return 1, 2, 4, or 8, matching the C++ scalar type `dt` represents.
 */
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

/**
 * @brief numpy dtype string (kind + itemsize) for a `DType`, e.g. `"f8"`, `"i4"`.
 * @param dt The dtype to convert.
 * @return A numpy-style struct format code understood by `numpy.dtype(...)`.
 */
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

/**
 * @brief A minimal typed, n-dimensional, row-major contiguous array.
 *
 * `NDArray` is either *owning* (holds its own `ByteBuf`, freed on
 * destruction) or a non-owning *view* over externally-managed memory
 * (`view_ != nullptr`); `is_view()` distinguishes the two, and `data()`
 * transparently returns whichever buffer is active. Views exist so the
 * write path can wrap a numpy array's memory directly (see
 * `bindings/np_conversions.hpp`'s `py_to_mesh`) without copying it into a
 * C++-owned buffer; `make_owned()` is the escape hatch for turning a view
 * into an owning copy when a buffer must outlive the memory it points to.
 * There is no reference counting: a view's caller is responsible for
 * keeping the underlying memory alive for the `NDArray`'s lifetime.
 */
class NDArray {
public:
    NDArray() = default;

    /**
     * @brief Constructs an owning array with a zero-initialized buffer.
     * @param dt Element dtype.
     * @param shape Row-major dimensions; total element count is their product.
     */
    NDArray(DType dt, std::vector<std::size_t> shape)
        : dtype_(dt), shape_(std::move(shape)) {
        const std::size_t nb = nbytes();
        owned_.resize(nb);  // uninitialised (no_init_allocator)
        std::memset(owned_.data(), 0, nb);  // explicit zero-fill
    }

    /**
     * @brief Constructs an owning array whose buffer is left *uninitialized*.
     *
     * Only safe for callers that immediately overwrite every byte — typical
     * uses are reader outputs (the whole buffer is about to be filled from
     * the parsed file) and cell-block reconstruction (e.g.
     * `detail::reconstruct_cells` in `vtk_cells.hpp`). Skips both the extra
     * allocator zero-fill and, more importantly, the cold first-touch page
     * faults a `memset` would otherwise incur on a fresh large allocation —
     * the same optimization numpy applies to its own `calloc`-avoidance path.
     * Prefer the two-argument constructor whenever the buffer might not be
     * fully overwritten.
     *
     * @param dt Element dtype.
     * @param shape Row-major dimensions; total element count is their product.
     * @return A new owning, uninitialized `NDArray`.
     */
    static NDArray uninit(DType dt, std::vector<std::size_t> shape) {
        NDArray a;
        a.dtype_ = dt;
        a.shape_ = std::move(shape);
        a.owned_.resize(a.nbytes());  // no memset
        return a;
    }

    /**
     * @brief Constructs a non-owning view over externally-owned row-major memory.
     *
     * Used to wrap a numpy array's buffer directly at the write boundary
     * (zero-copy): the C++ writer reads through `ptr` but never frees it.
     * @param dt Element dtype of the memory at `ptr`.
     * @param shape Row-major dimensions describing how to interpret `ptr`.
     * @param ptr Pointer to caller-owned memory; the caller must keep it
     *            alive for at least the lifetime of the returned `NDArray`
     *            (and of any `NDArray` copies/moves derived from it that
     *            remain a view).
     * @return A new non-owning `NDArray` view.
     */
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
    /** @brief Whether this array is a non-owning view (vs. owning its buffer). */
    bool is_view() const { return view_ != nullptr; }

    /** @brief Total element count (product of `shape()`), or 0 if `shape()` is empty. */
    std::size_t size() const {
        if (shape_.empty()) return 0;
        return std::accumulate(shape_.begin(), shape_.end(), std::size_t{1},
                               std::multiplies<std::size_t>());
    }
    /** @brief Total buffer size in bytes: `size() * dtype_size(dtype())`. */
    std::size_t nbytes() const { return size() * dtype_size(dtype_); }

    /** @brief Raw pointer to the active buffer (owned or view), for writing. */
    std::byte* data() { return view_ ? view_ : owned_.data(); }
    /** @brief Raw pointer to the active buffer (owned or view), read-only. */
    const std::byte* data() const { return view_ ? view_ : owned_.data(); }

    /**
     * @brief Changes the logical shape in place without touching the buffer.
     *
     * A no-op if the new shape's element count doesn't match the current
     * one (the mismatched reshape is silently ignored rather than throwing).
     * @param new_shape The desired row-major dimensions.
     */
    void reshape(std::vector<std::size_t> new_shape) {
        std::size_t n = new_shape.empty()
                            ? 0
                            : std::accumulate(new_shape.begin(), new_shape.end(),
                                              std::size_t{1}, std::multiplies<std::size_t>());
        if (n != size()) return;  // ignore inconsistent reshape
        shape_ = std::move(new_shape);
    }

    /**
     * @brief Turns a view into an owning copy in place; a no-op if already owning.
     *
     * Copies the viewed memory into a freshly-allocated owned buffer and
     * clears the view pointer. Used before handing a buffer's lifetime over
     * to Python via a capsule (`mesh_to_py`), where the destination `NDArray`
     * must actually own the memory it hands off.
     */
    void make_owned() {
        if (view_ == nullptr) return;
        const std::size_t nb = nbytes();
        ByteBuf buf;
        buf.resize(nb);  // uninitialised; fully overwritten by the memcpy below
        std::memcpy(buf.data(), view_, nb);
        owned_ = std::move(buf);
        view_ = nullptr;
    }

    /**
     * @brief Reinterprets the raw buffer as a `T*`. No dtype check is performed
     * — the caller must ensure `T` matches `dtype()`.
     * @tparam T The scalar type to view the buffer as.
     * @return Pointer to the first element, typed as `T`.
     */
    template <typename T>
    T* as() { return reinterpret_cast<T*>(data()); }
    /** @brief `const` overload of `as()`. */
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
