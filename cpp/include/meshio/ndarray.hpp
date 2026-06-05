#pragma once
//
// NDArray: a minimal typed, n-dimensional, contiguous (row-major) array.
//
// It is the storage primitive used by meshio::Mesh for points, cell
// connectivity, and all data fields. An NDArray either *owns* its buffer
// (the common case, e.g. data produced by a reader) or is a non-owning
// *view* over external memory (used to wrap numpy buffers zero-copy when
// writing). The binding layer converts between NDArray and numpy.

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <numeric>
#include <vector>

namespace meshio {

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
        owned_.resize(nbytes());
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

    // Turn a view into an owning copy (no-op if already owning). Used before
    // handing a buffer's lifetime to Python via a capsule.
    void make_owned() {
        if (view_ == nullptr) return;
        std::vector<std::byte> buf(nbytes());
        std::memcpy(buf.data(), view_, nbytes());
        owned_ = std::move(buf);
        view_ = nullptr;
    }

    template <typename T>
    T* as() { return reinterpret_cast<T*>(data()); }
    template <typename T>
    const T* as() const { return reinterpret_cast<const T*>(data()); }

private:
    DType dtype_ = DType::Float64;
    std::vector<std::size_t> shape_;
    std::vector<std::byte> owned_;
    std::byte* view_ = nullptr;
};

}  // namespace meshio
