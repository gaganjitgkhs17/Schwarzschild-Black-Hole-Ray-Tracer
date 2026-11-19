"""Version 0.10S — independent NumPy-only snapshot engine.

Equations and event logic copied from the archived 0.9A/0.9B source snapshot,
not from the user's potentially newer locally corrected 0.9A. Check with
physics_audit_v010s.py before adopting it. The radial camera ray remains
excluded as in the source. No SciPy, Matplotlib, or legacy-module imports.
"""
import numpy as np
M = 1.0
R_OBSERVER = 100.0
FOCAL_LENGTH = 100.0
HALF_WIDTH = 17.0
INCLINATION_DEG = 70.0
DISK_INNER = 6.0
DISK_OUTER = 20.0
STEP_SIZE = 0.008
MAX_STEPS = 6000
ROOT_ITERATIONS = 32
EXPOSURE = 2.5


def rhs(u, v):
    """
    Schwarzschild photon orbit equation:

        d²u/dphi² = 3*M*u² - u

    where:
        u = 1/r
        v = du/dphi
    """

    first_derivative = v

    second_derivative = (
        3.0 * M * u * u - u
    )

    return first_derivative, second_derivative


def rk4(u, v, h):
    """
    Advance photon trajectories using RK4.

    Works with individual values or NumPy arrays.
    """

    a1, b1 = rhs(u, v)

    a2, b2 = rhs(
        u + h * a1 / 2,
        v + h * b1 / 2
    )

    a3, b3 = rhs(
        u + h * a2 / 2,
        v + h * b2 / 2
    )

    a4, b4 = rhs(
        u + h * a3,
        v + h * b3
    )

    new_u = u + h * (
        a1 + 2*a2 + 2*a3 + a4
    ) / 6

    new_v = v + h * (
        b1 + 2*b2 + 2*b3 + b4
    ) / 6

    return new_u, new_v


def disk_plane_value(phi, az_z, inclination):
    """
    Evaluate the signed disk-plane condition.

    The disk normal is:
        n = (cos(i), 0, sin(i))

    A crossing occurs when the value changes sign.
    """

    return (
        np.cos(inclination) * np.cos(phi)
        +
        np.sin(inclination)
        * az_z
        * np.sin(phi)
    )


def refine_crossing(
    phi0,
    u0,
    v0,
    h,
    az_z,
    inclination
):
    """
    Refine a bracketed disk-plane crossing.

    Uses bisection to locate the crossing angle,
    followed by a fractional RK4 step to obtain
    the reciprocal radius at the crossing.
    """

    low = np.zeros_like(u0)

    high = np.full_like(
        u0,
        h
    )

    initial_sign = disk_plane_value(
        phi0,
        az_z,
        inclination
    )

    # Refine crossing location
    for iteration in range(ROOT_ITERATIONS):

        middle = (
            low + high
        ) / 2

        value = disk_plane_value(
            phi0 + middle,
            az_z,
            inclination
        )

        same_side = (
            initial_sign * value > 0
        )

        low = np.where(
            same_side,
            middle,
            low
        )

        high = np.where(
            same_side,
            high,
            middle
        )

    # Best estimate of angular crossing location
    crossing_step = (
        low + high
    ) / 2

    # Integrate the photon to this location
    u_hit, _ = rk4(
        u0,
        v0,
        crossing_step
    )

    phi_hit = (
        phi0 + crossing_step
    )

    return phi_hit, u_hit


def trace_disk(
    camera_x,
    camera_y,
    step_size=STEP_SIZE
):
    """
    Trace photons backward from the camera.

    Detect the first intersection with the disk
    between DISK_INNER and DISK_OUTER.

    Returns:
        hit_radius
        hit_angle
        steps_taken

    Pixels without disk intersections return NaN.
    """

    x = np.asarray(
        camera_x,
        dtype=float
    )

    y = np.asarray(
        camera_y,
        dtype=float
    )

    if x.shape != y.shape:
        raise ValueError(
            "Camera coordinate arrays must match."
        )

    if not np.isfinite(step_size) or step_size <= 0:
        raise ValueError("step_size must be positive and finite")
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        raise ValueError("camera coordinates must be finite")

    original_shape = x.shape

    # Flatten arrays for vectorized integration
    xx = x.ravel()
    yy = y.ravel()

    rho = np.hypot(
        xx,
        yy
    )

    # Metric factor at observer
    f_obs = (
        1.0 - 2.0 * M / R_OBSERVER
    )

    # Local camera viewing angle
    delta = np.arctan2(
        rho,
        FOCAL_LENGTH
    )

    # Photon impact parameters
    b = (
        R_OBSERVER
        * np.sin(delta)
        / np.sqrt(f_obs)
    )

    # Ray orientation in the camera plane
    az_z = np.divide(
        yy,
        rho,
        out=np.zeros_like(yy),
        where=rho > 0
    )

    inclination = np.deg2rad(
        INCLINATION_DEG
    )

    # Exclude the central radial ray
    active = b > 1e-12

    # Initial reciprocal radius
    u = np.full(
        xx.shape,
        1.0 / R_OBSERVER
    )

    # Initial radial derivative
    v = np.zeros_like(u)

    v[active] = np.sqrt(
        np.maximum(
            0.0,
            1.0 / b[active]**2
            - u[active]**2
            + 2.0 * M * u[active]**3
        )
    )

    phi = np.zeros_like(u)

    # Store intersection results
    hit_radius = np.full(
        xx.shape,
        np.nan
    )

    hit_angle = np.full(
        xx.shape,
        np.nan
    )

    # Numerical boundaries
    horizon_u = 1.0 / (2.0001 * M)
    sky_u = 1.0 / (400.0 * M)

    steps_taken = 0

    # ----------------------------------------------
    # INTEGRATE ALL ACTIVE PHOTONS
    # ----------------------------------------------

    for step in range(MAX_STEPS):

        indices = np.flatnonzero(active)

        if indices.size == 0:
            break

        steps_taken = step + 1

        # Current photon states
        u0 = u[indices]
        v0 = v[indices]
        p0 = phi[indices]

        # Advance trajectories
        u1, v1 = rk4(
            u0,
            v0,
            step_size
        )

        p1 = p0 + step_size

        # Evaluate disk plane before and after step
        z0 = disk_plane_value(
            p0,
            az_z[indices],
            inclination
        )

        z1 = disk_plane_value(
            p1,
            az_z[indices],
            inclination
        )

        # Find possible plane crossings
        crossing = (
            z0 * z1 <= 0
        )

        # ------------------------------------------
        # REFINE INTERSECTION LOCATIONS
        # ------------------------------------------

        if np.any(crossing):

            candidates = indices[crossing]

            p_hit, u_hit = refine_crossing(
                p0[crossing],
                u0[crossing],
                v0[crossing],
                step_size,
                az_z[candidates],
                inclination
            )

            # Convert reciprocal radius to radius
            radii = np.divide(
                1.0,
                u_hit,
                out=np.full_like(
                    u_hit,
                    np.inf
                ),
                where=u_hit > 0
            )

            # Check whether crossing is on disk
            valid = (
                (radii >= DISK_INNER)
                &
                (radii <= DISK_OUTER)
                &
                (u_hit < 1.0 / (2.0 * M))
            )

            selected = candidates[valid]

            # Record first disk intersection
            hit_radius[selected] = radii[valid]
            hit_angle[selected] = p_hit[valid]

            # Stop rays after disk hit
            active[selected] = False

        # ------------------------------------------
        # UPDATE REMAINING PHOTONS
        # ------------------------------------------

        remaining = active[indices]

        remaining_indices = indices[remaining]

        u[remaining_indices] = u1[remaining]
        v[remaining_indices] = v1[remaining]
        phi[remaining_indices] = p1[remaining]

        # Stop rays at horizon or distant sky
        done = (
            (u[remaining_indices] >= horizon_u)
            |
            (
                (u[remaining_indices] <= sky_u)
                &
                (v[remaining_indices] < 0)
            )
            |
            (~np.isfinite(u[remaining_indices]))
        )

        active[
            remaining_indices[done]
        ] = False

    if np.any(active):
        print(
            f"Warning: {np.count_nonzero(active)} "
            f"rays reached MAX_STEPS."
        )

    return (
        hit_radius.reshape(original_shape),
        hit_angle.reshape(original_shape),
        steps_taken
    )


def calculate_impact_parameter(x, y):
    """
    Convert local camera-plane coordinates to
    Schwarzschild photon impact parameter.

    b = R*sin(delta)/sqrt(1 - 2M/R)
    """

    rho = np.hypot(x, y)

    delta = np.arctan2(
        rho,
        FOCAL_LENGTH
    )

    f_obs = 1.0 - 2.0 * M / R_OBSERVER

    b = (
        R_OBSERVER
        * np.sin(delta)
        / np.sqrt(f_obs)
    )

    return b


def photon_angular_momentum(x, y):
    """
    Calculate the signed photon angular momentum
    about the disk's rotation axis.

    The camera is on the positive x-axis.

    Camera horizontal coordinate maps to the
    positive y-direction for backward-traced rays.

    Disk axis:
        n = (cos(i), 0, sin(i))

    The physical photon travels from the disk
    toward the camera, opposite to the backward
    ray-tracing direction.

    Therefore:
        lambda = -b*sin(i)*(camera_x/rho)

    This sign sets the disk's rotation convention.
    """

    b = calculate_impact_parameter(x, y)

    rho = np.hypot(x, y)

    inclination = np.deg2rad(
        INCLINATION_DEG
    )

    direction_y = np.divide(
        x,
        rho,
        out=np.zeros_like(rho, dtype=float),
        where=rho > 0
    )

    angular_momentum = (
        -b
        * np.sin(inclination)
        * direction_y
    )

    return angular_momentum


def orbital_angular_velocity(radius):
    """
    Circular geodesic angular velocity:

        Omega = sqrt(M/r^3)

    The disk begins at r = 6M, the
    Schwarzschild innermost stable circular orbit.
    """

    return np.sqrt(
        M / radius**3
    )


def frequency_shift(radius, angular_momentum):
    """
    Calculate observed-to-emitted frequency ratio.

    g = sqrt(1 - 3M/r) /
        [sqrt(1 - 2M/R_obs) * (1 - Omega*lambda)]

    The observer is static at R_OBSERVER.
    The emitting disk follows circular geodesics.
    """

    omega = orbital_angular_velocity(
        radius
    )

    f_observer = (
        1.0 - 2.0 * M / R_OBSERVER
    )

    emitter_factor = np.sqrt(
        1.0 - 3.0 * M / radius
    )

    doppler_denominator = (
        1.0 - omega * angular_momentum
    )

    if np.any(doppler_denominator <= 0):
        raise ValueError(
            "Invalid Doppler denominator. "
            "Check photon angular momentum."
        )

    g = emitter_factor / (
        np.sqrt(f_observer)
        * doppler_denominator
    )

    return g


def emitted_intensity(radius):
    """
    Simplified radial emission profile.

    Inner regions emit more strongly than
    outer regions.

    This is an illustrative emission law,
    not a complete accretion-disk model.
    """

    return (
        DISK_INNER / radius
    ) ** 3


def emitted_colour(radius):
    """
    Assign illustrative linear RGB colour
    to emitted disk radiation.

    Inner disk is warmer-white.
    Outer disk is red-orange.

    This colour assignment is not a
    calibrated thermal spectrum.
    """

    t = np.clip(
        (radius - DISK_INNER)
        / (DISK_OUTER - DISK_INNER),
        0.0,
        1.0
    )

    red = np.ones_like(t)

    green = (
        0.38 + 0.55 * (1.0 - t)
    )

    blue = (
        0.07 + 0.65 * (1.0 - t)**2
    )

    return np.stack(
        (red, green, blue),
        axis=-1
    )


def tone_map(linear_rgb):
    """
    Convert nonnegative linear RGB intensities
    into display values between 0 and 1.

    This prevents abrupt clipping of bright areas.
    """

    return 1.0 - np.exp(
        -EXPOSURE * np.maximum(linear_rgb, 0.0)
    )



def trace_radiance(x, y, step_size=STEP_SIZE):
    """Return linear RGB and disk-hit flags; no tone mapping."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.shape != y.shape: raise ValueError('Coordinate shape mismatch')
    radii, _, _ = trace_disk(x, y, step_size=step_size)
    hits = np.isfinite(radii)
    linear = np.zeros(x.shape + (3,), dtype=np.float64)
    if np.any(hits):
        lam = photon_angular_momentum(x, y)
        shift = frequency_shift(radii[hits], lam[hits])
        linear[hits] = (emitted_intensity(radii[hits]) * shift**4)[:, None] * emitted_colour(radii[hits])
    if not np.isfinite(linear).all() or (linear < 0).any():
        raise RuntimeError('Invalid ray intensity')
    return linear, hits


def render_linear(size=300, samples=4, tile_rows=8, step_size=STEP_SIZE, verbose=True):
    """Return (linear RGB, per-pixel hit fraction), tiled, exactly sampled."""
    if not all(isinstance(n, int) and not isinstance(n, bool) and n > 0 for n in (size,samples,tile_rows)):
        raise ValueError('size/samples/tile_rows must be positive integers')
    pitch = 2 * HALF_WIDTH / size
    xs = -HALF_WIDTH + (np.arange(size*samples,dtype=float)+.5)*pitch/samples
    rgb = np.empty((size,size,3),dtype=float)
    coverage = np.empty((size,size),dtype=float)
    for row0 in range(0,size,tile_rows):
        row1 = min(size,row0+tile_rows)
        ys = -HALF_WIDTH + (np.arange(row0*samples,row1*samples,dtype=float)+.5)*pitch/samples
        xx,yy = np.meshgrid(xs,ys)
        vals,hits=trace_radiance(xx,yy,step_size)
        rgb[row0:row1] = vals.reshape(row1-row0,samples,size,samples,3).mean(axis=(1,3))
        coverage[row0:row1] = hits.reshape(row1-row0,samples,size,samples).mean(axis=(1,3))
        if verbose: print(f'{row1:4d}/{size} rows',flush=True)
    return rgb, coverage


def camera_brightness_integral(linear):
    return float(np.mean(linear, axis=2).sum() * (2*HALF_WIDTH/linear.shape[0])**2)
