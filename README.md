# ⚡ Scalable Cloud-in-Cell Interpolation (OpenMP-C)

Efficient parallel implementation of **bilinear interpolation** from scattered 2D points onto a structured mesh using **OpenMP**.This project demonstrates parallel computing techniques, cache-aware optimization, and performance profiling on modern multi-core systems.

---

## 📊 Problem Statement

Given a set of scattered 2D points $(x_i, y_i)$, interpolate values onto a structured $X_CELLS \times Y_CELLS$ mesh using the **Cloud-in-Cell (CIC)** method — a bilinear interpolation approach widely used in simulations, graphics, and numerical modeling.

---

## 🚀 Key Highlights

* ✅ Optimized OpenMP parallelization
* ↻ Multiple iterations with timing support
* 🧮 Guided scheduling for load balancing
* ❌ Thread-safe local grid buffers (no race conditions)

---

## ⚙️ How It Works

### 📅 Input

A binary file containing:

* Grid dimensions: `X_CELLS`, `Y_CELLS`
* Number of particles
* Number of iterations
* Double-precision coordinates $(x, y)$ of each particle

### 🔧 Components

| Function             | Description                                                  |
| -------------------- | ------------------------------------------------------------ |
| `main()`             | Reads input, controls timing and iterations                  |
| `loadParticleData()` | Efficient bulk read of particle positions with parallel copy |
| `distributeCharge()` | CIC interpolation using private buffers & guided scheduling  |
| `dumpGridToFile()`   | Writes interpolated mesh to `Mesh1.out` in readable format   |

---

## 🔋 Parallelization Strategy

The `distributeCharge()` function employs multiple strategies to ensure efficient and scalable parallelization:

* **Thread-local grid buffers**: Each thread maintains its own private buffer to accumulate interpolated values. This prevents race conditions and removes the need for atomic operations or locks during grid updates.

* **Guided scheduling for particle loop**: The particle processing loop is parallelized using `schedule(guided, 256)`. This means larger chunks of particles are assigned at first, followed by smaller chunks. It improves load balancing, especially when particle distribution is uneven or access costs vary.

* **Static scheduling for initialization and reduction**: For loops like grid initialization and final reduction (summing all thread-local grids into the global grid), `schedule(static)` is used. Static scheduling divides work evenly among threads and is optimal when each iteration has similar cost, reducing scheduling overhead.

* **Memory locality**: Since each thread accesses its own buffer, memory accesses are contiguous and cache-friendly. This improves spatial locality and reduces cache contention and false sharing.

* **Final reduction phase**: After interpolation, a parallel reduction phase efficiently combines thread-local results into the final global grid. This step is parallelized too, avoiding bottlenecks.

Overall, the combination of thread-private buffers, guided scheduling for dynamic work distribution, and static scheduling for uniform tasks leads to a fast, scalable, and race-free implementation of the Cloud-in-Cell interpolation.

---

## 📊 Performance Comparison

### 🔢 Speedup Table: All Cases vs Number of Threads

| Threads | Case 1 (0.9M) | Case 2 (5M) | Case 3 (3.6M) | Case 4 (20M) | Case 5 (14M) |
| ------- | ------------- | ----------- | ------------- | ------------ | ------------ |
| 1       | 1.0           | 1.0         | 1.0           | 1.0          | 1.0          |
| 2       | 2.0           | 2.1         | 1.5           | 1.4          | 1.3          |
| 3       | 2.3           | 2.8         | 2.0           | 2.1          | 1.9          |
| 4       | 2.7           | 3.7         | 2.2           | 2.6          | 2.3          |
| 5       | 3.2           | 4.3         | 2.6           | 2.9          | 2.6          |
| 6       | 4.3           | 5.1         | 3.0           | 3.7          | 3.1          |
| 7       | 5.1           | 5.8         | 3.3           | 3.9          | 3.5          |
| 8       | 5.8           | 6.5         | 3.6           | 4.7          | 4.0          |
| 9       | 6.3           | 7.0         | 3.8           | 5.2          | 4.3          |
| 10      | 6.9           | 7.7         | 4.1           | 5.7          | 4.4          |
| 11      | 7.4           | 8.2         | 4.4           | 5.7          | 4.5          |
| 12      | 7.9           | 8.5         | 4.7           | 6.3          | 4.6          |

---

## 📅 Build & Run

### ✅ Compile

```bash
gcc -fopenmp src/cic_parallel_optimized.c -o cic_parallel_optimized
```

### ▶️ Run

```bash
./cic_parallel_optimized  sample_input_file.bin  8
```

* sample_input_file.bin is the input file generated using src/input_fileMaker.c
* `8` is the number of OpenMP threads
* Output is written to `Mesh1.out`

---

## 📜 Report

Detailed analysis, hardware specifications, graphs, and theoretical insight are available in:
[`Implementation_report.pdf`](./Implementation_report.pdf)

---
