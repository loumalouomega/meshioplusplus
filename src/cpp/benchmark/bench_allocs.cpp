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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
// System includes
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <new>

// Project includes
#include "bench_allocs.hpp"

namespace {

struct BenchAllocHeader {
    void* mpBase;
    std::size_t mSize;
    std::uint64_t mEpoch;
};

std::atomic<bool> g_alloc_active{false};
std::atomic<std::uint64_t> g_alloc_epoch{1};
std::atomic<std::uint64_t> g_alloc_count{0};
std::atomic<std::uint64_t> g_alloc_total_bytes{0};
std::atomic<std::uint64_t> g_alloc_live_bytes{0};
std::atomic<std::uint64_t> g_alloc_peak_bytes{0};

void* bench_alloc_raw(std::size_t Size, std::size_t Align) noexcept {
    std::size_t align = Align < alignof(std::max_align_t) ? alignof(std::max_align_t) : Align;
    if (align < alignof(BenchAllocHeader))
        align = alignof(BenchAllocHeader);
    if (Size > std::numeric_limits<std::size_t>::max() - align - sizeof(BenchAllocHeader))
        return nullptr;
    const std::size_t total = Size + align + sizeof(BenchAllocHeader);
    void* raw = std::malloc(total);
    if (!raw)
        return nullptr;
    const auto base = reinterpret_cast<std::uintptr_t>(raw) + sizeof(BenchAllocHeader);
    const auto aligned = (base + (align - 1)) & ~(align - 1);
    auto* hdr = reinterpret_cast<BenchAllocHeader*>(aligned) - 1;
    hdr->mpBase = raw;
    hdr->mSize = Size;
    if (g_alloc_active.load(std::memory_order_acquire)) {
        hdr->mEpoch = g_alloc_epoch.load(std::memory_order_relaxed);
        g_alloc_count.fetch_add(1, std::memory_order_relaxed);
        g_alloc_total_bytes.fetch_add(Size, std::memory_order_relaxed);
        const std::uint64_t live =
            g_alloc_live_bytes.fetch_add(Size, std::memory_order_relaxed) + Size;
        std::uint64_t peak = g_alloc_peak_bytes.load(std::memory_order_relaxed);
        while (live > peak &&
               !g_alloc_peak_bytes.compare_exchange_weak(peak, live, std::memory_order_relaxed)) {
        }
    } else {
        hdr->mEpoch = 0;
    }
    return reinterpret_cast<void*>(aligned);
}

void* bench_alloc_throwing(std::size_t Size, std::size_t Align) {
    while (true) {
        if (void* p = bench_alloc_raw(Size, Align))
            return p;
        std::new_handler handler = std::get_new_handler();
        if (!handler)
            throw std::bad_alloc();
        handler();
    }
}

void bench_free_raw(void* pPtr) noexcept {
    if (!pPtr)
        return;
    const BenchAllocHeader hdr = *(reinterpret_cast<const BenchAllocHeader*>(pPtr) - 1);
    if (hdr.mEpoch != 0 && g_alloc_active.load(std::memory_order_acquire) &&
        hdr.mEpoch == g_alloc_epoch.load(std::memory_order_relaxed)) {
        g_alloc_live_bytes.fetch_sub(hdr.mSize, std::memory_order_relaxed);
    }
    std::free(hdr.mpBase);
}

}  // namespace

void* operator new(std::size_t Size) {
    return bench_alloc_throwing(Size, alignof(std::max_align_t));
}

void* operator new[](std::size_t Size) {
    return bench_alloc_throwing(Size, alignof(std::max_align_t));
}

void* operator new(std::size_t Size, const std::nothrow_t&) noexcept {
    return bench_alloc_raw(Size, alignof(std::max_align_t));
}

void* operator new[](std::size_t Size, const std::nothrow_t&) noexcept {
    return bench_alloc_raw(Size, alignof(std::max_align_t));
}

void* operator new(std::size_t Size, std::align_val_t Align) {
    return bench_alloc_throwing(Size, static_cast<std::size_t>(Align));
}

void* operator new[](std::size_t Size, std::align_val_t Align) {
    return bench_alloc_throwing(Size, static_cast<std::size_t>(Align));
}

void* operator new(std::size_t Size, std::align_val_t Align, const std::nothrow_t&) noexcept {
    return bench_alloc_raw(Size, static_cast<std::size_t>(Align));
}

void* operator new[](std::size_t Size, std::align_val_t Align, const std::nothrow_t&) noexcept {
    return bench_alloc_raw(Size, static_cast<std::size_t>(Align));
}

void operator delete(void* pPtr) noexcept {
    bench_free_raw(pPtr);
}
void operator delete[](void* pPtr) noexcept {
    bench_free_raw(pPtr);
}
void operator delete(void* pPtr, std::size_t) noexcept {
    bench_free_raw(pPtr);
}
void operator delete[](void* pPtr, std::size_t) noexcept {
    bench_free_raw(pPtr);
}
void operator delete(void* pPtr, const std::nothrow_t&) noexcept {
    bench_free_raw(pPtr);
}
void operator delete[](void* pPtr, const std::nothrow_t&) noexcept {
    bench_free_raw(pPtr);
}
void operator delete(void* pPtr, std::align_val_t) noexcept {
    bench_free_raw(pPtr);
}
void operator delete[](void* pPtr, std::align_val_t) noexcept {
    bench_free_raw(pPtr);
}
void operator delete(void* pPtr, std::size_t, std::align_val_t) noexcept {
    bench_free_raw(pPtr);
}
void operator delete[](void* pPtr, std::size_t, std::align_val_t) noexcept {
    bench_free_raw(pPtr);
}
void operator delete(void* pPtr, std::align_val_t, const std::nothrow_t&) noexcept {
    bench_free_raw(pPtr);
}
void operator delete[](void* pPtr, std::align_val_t, const std::nothrow_t&) noexcept {
    bench_free_raw(pPtr);
}

namespace meshioplusplus::bench {

void alloc_counter_start() noexcept {
    g_alloc_count.store(0, std::memory_order_relaxed);
    g_alloc_total_bytes.store(0, std::memory_order_relaxed);
    g_alloc_live_bytes.store(0, std::memory_order_relaxed);
    g_alloc_peak_bytes.store(0, std::memory_order_relaxed);
    g_alloc_epoch.fetch_add(1, std::memory_order_relaxed);
    g_alloc_active.store(true, std::memory_order_release);
}

AllocStats alloc_counter_stop() noexcept {
    g_alloc_active.store(false, std::memory_order_release);
    return {g_alloc_count.load(std::memory_order_relaxed),
            g_alloc_total_bytes.load(std::memory_order_relaxed),
            g_alloc_peak_bytes.load(std::memory_order_relaxed)};
}

}  // namespace meshioplusplus::bench
