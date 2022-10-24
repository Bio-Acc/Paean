#ifndef PAEAN_BISECT_H
#define PAEAN_BISECT_H

#include <stdint.h>

int bisect_left(uint64_t *array, uint32_t numOfEntry, uint64_t x) {
    uint32_t left = 0, right = numOfEntry, mid;
    while (left < right) {
        mid = (left + right) / 2;
        if (array[mid] < x)
            left = mid + 1;
        else
            right = mid;
    }
    return left;
}

int bisect_right(uint64_t *array, uint32_t numOfEntry, uint64_t x) {
    uint32_t left = 0, right = numOfEntry, mid;
    while (left < right) {
        mid = (left + right) / 2;
        if (x < array[mid])
            right = mid;
        else
            left = mid + 1;
    }
    return left;
}
int find_coincide(uint64_t key, uint64_t *entries, uint32_t numOfEntry) {
    // find leftmost entry exactly equal to key
    uint32_t i = bisect_left(entries, numOfEntry, key);
    if (i != numOfEntry && entries[i] == key)
        return i;
    else
        return -1;
}

void find_inner_contain(uint64_t k_start, uint64_t k_end,
                        uint64_t *d_starts,
                        uint32_t numOfEntry,
                        uint32_t *d_assist_starts,
                        uint32_t *d_assist_ends, uint32_t id) {
    int i;
    // find leftmost entry whose start >= k_start
    i = bisect_left(d_starts, numOfEntry, k_start);
    if (i != numOfEntry)
        d_assist_starts[id] = i;
    else {
        d_assist_starts[id] = d_assist_ends[id] = 0;
        return;
    }
    // find leftmost entry whose start >= k_end
    i = bisect_left(d_starts, numOfEntry, k_end);
    if (i) {
        d_assist_ends[id] = i;
    } else {
        d_assist_starts[id] = d_assist_ends[id] = 0;
    }
}

int count_inner_contain(uint64_t k_start, uint64_t k_end,
                        uint64_t *d_starts, uint64_t *d_ends,
                        uint32_t numOfEntry) {
    int i, left;
    // find leftmost entry whose start >= k_start
    i = bisect_left(d_starts, numOfEntry, k_start);
    if (i != numOfEntry)
        left = i;
    else {
        return 0;
    }
    // find rightmost entry whose start < k_end
    i = bisect_left(d_starts, numOfEntry, k_end);
    if (i) {
        int c = 0;
        while (--i >= left) {
            if (d_ends[i] <= k_end)
                c++;
        }
        return c;
    } else {
        return 0;
    }
}

int count_outer_contain(uint64_t k_start, uint64_t k_end,
                        uint64_t *d_starts, uint64_t *d_ends,
                        uint32_t numOfEntry) {
    int left, right;
    // find leftmost entry whose start >= k_start
    left = bisect_left(d_starts, numOfEntry, k_start);
    // find rightmost entry whose end <= k_end
    right = bisect_right(d_ends, numOfEntry, k_end);
    return (left - right);
}

int count_overlap(uint64_t k_start, uint64_t k_end,
                  uint64_t *d_starts, uint64_t *d_ends,
                  uint32_t numOfEntry) {
    int left, right;
    // find leftmost entry whose start > k_end
    left = bisect_right(d_starts, numOfEntry, k_end);
    // // find rightmost entry whose end < k_start
    right = bisect_left(d_ends, numOfEntry, k_start);
    return (left - right);
}

#endif // PAEAN_BISECT_H