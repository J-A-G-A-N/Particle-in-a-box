import time
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from matplotlib.animation import FFMpegFileWriter
from numba import njit, prange


def timer_dec(base_func):
    def extended_func(*args, **kwargs):
        start = time.time()
        result = base_func(*args, **kwargs)
        end = time.time()
        enlapsed_s = end - start
        print(f"enlapsed_us of {base_func.__name__}: {enlapsed_s:.2f}s ")
        return result

    return extended_func


def build_grid(N: int, L: float):
    x = np.linspace(-L, L, N)
    dx = x[1] - x[0]
    return x, dx


def finite_barrier(x, x1, x2, V0):
    V = np.zeros_like(x)
    mask = (x >= x1) & (x <= x2)
    V[mask] = V0
    return V


def build_hamiltonian(N, dx, potenetial):
    B = -1 / (2 * dx**2)
    A = (1 / dx**2) + potenetial

    main = np.diag(A)
    upper = np.diag(np.ones(N - 1) * B, 1)
    lower = np.diag(np.ones(N - 1) * B, -1)
    return main + upper + lower


@njit
def psi_t(energies, basis, coeffs, t):
    hbar = 1.0
    phase = np.exp(-1j * energies * t / hbar)
    return np.sum((coeffs * phase)[:, None] * basis, axis=0)


def solve_eigh(H):
    eig_val, eig_vec = sp.linalg.eigh(H)
    idx = np.argsort(eig_val)
    return eig_val[idx], eig_vec[:, idx].T


@timer_dec
@njit(parallel=True)
def compute_density_fast(time, energies, basis, coeffs):
    T = len(time)
    M = basis.shape[1]
    density = np.empty((T, M), dtype=np.float64)
    for i in prange(T):
        psi_val = psi_t(energies, basis, coeffs, time[i])
        density[i, :] = np.abs(psi_val) ** 2
    return density


@timer_dec
def compute_density(time, energies, basis, coeffs):
    return [np.abs(psi_t(energies, basis, coeffs, t)) ** 2 for t in time]


def normalize(psi, dx):
    return psi / np.sqrt(np.sum(np.abs(psi) ** 2) * dx)


def compute_coefficients(basis, psi0, x):
    return np.trapezoid(np.conj(basis) * psi0[None, :], x, axis=1)


def energies_in_ev(engeries, L):
    hbar = sp.constants.hbar
    m_e = sp.constants.m_e
    ev = sp.constants.e
    L_si = L * 1e-9
    E_scale_ev = (hbar**2 / (m_e * L_si**2)) / ev
    return engeries * E_scale_ev


@timer_dec
def solve_numerical(N, dx, n_modes, barrier_start, barrier_end, barrier_height, x):
    V = finite_barrier(x, barrier_start, barrier_end, barrier_height)
    H = build_hamiltonian(N, dx, V)
    eigen_values, eigen_vectors = solve_eigh(
        H
    )  # Sorted and Normalized eigen_values and eigen_vectors
    state_vec = np.array([normalize(vec, dx) for vec in eigen_vectors])
    return eigen_values[:n_modes], state_vec[:n_modes]


def initial_superposition(psi_modes, n_super, dx):
    psi0 = np.sum(psi_modes[:n_super], axis=0)
    return normalize(psi0, dx)


def plot_energy_vs_n(n_values, numerical_energies):
    if n_values > len(numerical_energies):
        raise ValueError("n_values exceeds available data")

    n = np.arange(1, n_values + 1)
    x = n

    num = numerical_energies[:n_values]

    plt.figure()

    plt.plot(x, num, "-", label="Numerical")

    plt.xlabel(r"$n$")
    plt.ylabel("Energy")
    plt.title("Energy vs $n$")
    plt.legend()
    plt.grid(True)
    plt.savefig("results/finite_barrier/energy_vs_n.png")


@timer_dec
def animate_density_video(x, density, outname, fps=30):
    ymax = np.max(density)
    fig, ax = plt.subplots()
    (line,) = ax.plot([], [], color="blue", label="Numerical")
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(0, ymax)
    ax.set_xlabel(r"$x / L$")
    ax.set_ylabel(r"$|\psi(x,t)|^2$")
    ax.legend(loc="upper right")
    writer = FFMpegFileWriter(fps, codec="gif")
    with writer.saving(fig, outname, dpi=200):
        for i, num in enumerate(density):
            line.set_data(x, num)
            ax.set_title(f"Time step:{i:04d}")
            writer.grab_frame()
    plt.close(fig)

def gaussian_packet(x, x0, sigma, k0):
    psi = (
        (1 / sigma * np.sqrt(2 * np.pi))
        * np.exp(-((x - x0) ** 2) / (2 * sigma**2))
        * np.exp(1j * k0 * x)
    )
    return psi


def main():
    L = 15
    N = 5000
    fps = 25
    n_modes = 500
    barrier_start = 0
    barrier_end = 0.2
    barrier_height = 5
    k0 = 5
    sigma = 0.8
    x0 = -10
    x, dx = build_grid(N, L)
    # -------------------------------------- Solve Finite Differentiation Method ------------------
    e, v = solve_numerical(
        N, dx, n_modes, barrier_start, barrier_end, barrier_height, x
    )

    psi0_numerical = normalize(gaussian_packet(x, x0, sigma, k0), dx)
    # -------------------------------------- Compute Co-efficients for Theoretical and Numerical---
    c_n = compute_coefficients(v, psi0_numerical, x)
    # -------------------------------------- Time Evolution ---------------------------------------
    frames = 30 * fps
    t_sim = 2 * L / k0
    time = np.linspace(0, t_sim * 4, frames)
    # density_numerical = compute_density(time, e, v, c_n)
    # animate_density_video(
    #     x, density_numerical, "results/finite_barrier/density.gif", fps
    # )
    plot_energy_vs_n(n_modes, e)
    print("Done")


if __name__ == "__main__":
    main()
