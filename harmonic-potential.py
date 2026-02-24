import os
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp


def build_grid(N, L):
    x = np.linspace(-L, L, N)
    dx = x[1] - x[0]
    return x, dx


def harmonic_potential(x, omega):
    return 0.5 * omega**2 * x**2


def build_hamiltonian(N, dx, potenetial):
    B = -1 / (2 * dx**2)
    A = (1 / dx**2) + potenetial

    main = np.diag(A)
    upper = np.diag(np.ones(N - 1) * B, 1)
    lower = np.diag(np.ones(N - 1) * B, -1)
    return main + upper + lower


def psi_t(eig_val, psi, c, t):
    hbar = 1
    phase = np.exp(-1j * eig_val * t / hbar)
    return np.sum((c * phase)[:, None] * psi, axis=0)


def solve_eigh(H):
    eig_val, eig_vec = sp.linalg.eigh(H)
    idx = np.argsort(eig_val)
    return eig_val[idx], eig_vec[:, idx].T


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


def solve_numerical(N, dx, n_modes, omega, x):
    V = harmonic_potential(x, omega)
    H = build_hamiltonian(N, dx, V)
    eigen_values, eigen_vectors = solve_eigh(
        H
    )  # Sorted and Normalized eigen_values and eigen_vectors
    state_vec = np.array([normalize(vec, dx) for vec in eigen_vectors])
    return eigen_values[:n_modes], state_vec[:n_modes]


def analytical_wavefunction(n_modes, omega, x):
    x = np.asarray(x)
    xi = np.sqrt(omega) * x
    pref = (omega / np.pi) ** 0.25

    states = np.empty((n_modes, x.size))

    exp_factor = np.exp(-0.5 * omega * x**2)

    for n in range(n_modes):
        norm = pref / np.sqrt(2**n * sp.special.factorial(n))
        states[n] = norm * sp.special.eval_hermite(n, xi) * exp_factor

    return states


def analytical_energies(n_modes, omega):
    n = np.arange(1, n_modes + 1)
    return omega * (n + 0.5)


def analytical_solution(n_modes, omega, x):
    energies = analytical_energies(n_modes, omega)
    states = analytical_wavefunction(n_modes, omega, x)
    print(states.shape)
    return energies, states


def fix_phase_shift(n_modes, psi_theory, psi_num, x):
    for i in range(n_modes):
        overlap = np.trapezoid(psi_num[i] * psi_theory[i], x)
        if overlap < 0:
            psi_num[i] *= -1
    return psi_num


def initial_superposition(psi_modes, n_super, dx):
    psi0 = np.sum(psi_modes[:n_super], axis=0)
    return normalize(psi0, dx)


def plot_energy_vs_n_sq(n_values, numerical_energies, theoretical_energies):
    if len(numerical_energies) != len(theoretical_energies):
        raise ValueError("Energy arrays must have same length")

    if n_values > len(numerical_energies):
        raise ValueError("n_values exceeds available data")

    n = np.arange(1, n_values + 1)
    x = n**2

    num = numerical_energies[:n_values]
    theo = theoretical_energies[:n_values]

    error = np.abs(num - theo)

    plt.figure()

    plt.errorbar(x, num, yerr=error, fmt="o", label="Numerical", capsize=4)

    plt.plot(x, theo, "-", label="Theoretical")

    plt.xlabel(r"$n^2$")
    plt.ylabel("Energy")
    plt.title("Energy vs $n^2$ with Error Bars")
    plt.legend()
    plt.grid(True)
    plt.savefig("results/harmonic/energy_vs_n^2.png")


def animate_density(x, density_theory, density_numerical, frames_outdir):
    ymax = max(np.max(density_theory), np.max(density_numerical))
    for i, (th, num) in enumerate(zip(density_theory, density_numerical)):
        plt.figure()

        plt.plot(x, th, label="Theory", color="red", linewidth=3)
        plt.plot(x, num, "--", label="Numerical", color="blue", linewidth=3)

        plt.xlim(x.min(), x.max())
        plt.ylim(0, ymax)

        plt.xlabel("x/L")
        plt.ylabel(r"$|\psi|^2$")
        plt.legend()
        plt.title(f"Frame {i:04d}")

        filename = os.path.join(frames_outdir, f"frame_{i:04d}.png")
        plt.savefig(filename, dpi=150)
        plt.close()

    print("Frames saved.")


def animate_density_numerical(x, density_numerical, frames_outdir):
    ymax = np.max(density_numerical)
    for i, num in enumerate(density_numerical):
        plt.figure()

        plt.plot(x, num, "--", label="Numerical", color="blue", linewidth=3)

        plt.xlim(x.min(), x.max())
        plt.ylim(0, ymax)

        plt.xlabel("x/L")
        plt.ylabel(r"$|\psi|^2$")
        plt.legend()
        plt.title(f"Frame {i:04d}")

        filename = os.path.join(frames_outdir, f"frame_{i:04d}.png")
        plt.savefig(filename, dpi=150)
        plt.close()

    print("Frames saved.")


def setup_results_dir():
    frames_dir = "results/harmonic/frames"
    os.makedirs("results/harmonic", exist_ok=True)
    os.makedirs(frames_dir, exist_ok=True)


def ffmpeg_make_video_and_gif(
    frames_pattern="results/harmonic/frames/frame_%04d.png",
    mp4_name="results/harmonic/density.mp4",
    gif_name="results/harmonic/density.gif",
    fps=30,
    scale_width=800,
):
    try:
        cmd_mp4 = [
            "ffmpeg",
            "-y",
            "-framerate",
            str(fps),
            "-i",
            frames_pattern,
            "-c:v",
            "libx264",
            "-pix_fmt",
            "yuv420p",
            mp4_name,
        ]

        subprocess.run(cmd_mp4, check=True)

        palette = "palette.png"
        cmd_palette = [
            "ffmpeg",
            "-y",
            "-i",
            mp4_name,
            "-vf",
            f"fps={fps},scale={scale_width}:-1:flags=lanczos,palettegen",
            palette,
        ]

        subprocess.run(cmd_palette, check=True)

        cmd_gif = [
            "ffmpeg",
            "-y",
            "-i",
            mp4_name,
            "-i",
            palette,
            "-lavfi",
            f"fps={fps},scale={scale_width}:-1:flags=lanczos[x];[x][1:v]paletteuse",
            gif_name,
        ]

        subprocess.run(cmd_gif, check=True)

        if os.path.exists(palette):
            os.remove(palette)

    except FileNotFoundError:
        print("FFmpeg not found. Install ffmpeg first.")
    except subprocess.CalledProcessError:
        print("FFmpeg command failed.")


def main():
    L = 5
    N = 4000
    fps = 25
    n_super = 2
    n_modes = 50
    omega = 50
    setup_results_dir()
    x, dx = build_grid(N, L)

    # -------------------------------------- Sovle Finite Differentiation Method ------------------
    e, v = solve_numerical(N, dx, n_modes, omega, x)

    # -------------------------------------- Calculated Theoritical values ------------------------
    ae, av = analytical_solution(n_modes, omega, x)
    # -------------------------------------- Account for phase Shift ------------------------------
    v = fix_phase_shift(n_modes, av, v, x)
    # -------------------------------------- Intial Wave function Theoritical and Numerical -------
    psi0_numerical = initial_superposition(v, n_super, dx)
    psi0_analytical = initial_superposition(av, n_super, dx)
    # -------------------------------------- Compute Co-efficients for Theoritical and Numerical---
    c_n = compute_coefficients(v, psi0_numerical, x)
    a_n = compute_coefficients(av, psi0_analytical, x)

    # -------------------------------------- Time Evolution ---------------------------------------
    frames = 5 * fps
    time = np.linspace(0, (2 * np.pi) / omega, frames)
    density_numerical = compute_density(time, e, v, c_n)
    density_theory = compute_density(time, ae, av, a_n)
    plot_energy_vs_n_sq(n_modes, e, ae)
    animate_density(x, density_numerical, density_theory, "results/harmonic/frames")
    ffmpeg_make_video_and_gif(fps=fps)
    print("Done")


if __name__ == "__main__":
    main()
