# Euler Parameter Demo
# Vibe coded via Gemini

import pyvista as pv
import numpy as np

# --- 1. Math Functions ---


def quaternion_to_dcm(q):
    """
    Converts a 4-element quaternion (q0, q1, q2, q3)
    into a 3x3 Direction Cosine Matrix (DCM).
    q0 is the scalar part.
    """
    q0, q1, q2, q3 = q

    R00 = q0**2 + q1**2 - q2**2 - q3**2
    R01 = 2 * (q1 * q2 - q0 * q3)
    R02 = 2 * (q1 * q3 + q0 * q2)
    R10 = 2 * (q1 * q2 + q0 * q3)
    R11 = q0**2 - q1**2 + q2**2 - q3**2
    R12 = 2 * (q2 * q3 - q0 * q1)
    R20 = 2 * (q1 * q3 - q0 * q2)
    R21 = 2 * (q2 * q3 + q0 * q1)
    R22 = q0**2 - q1**2 - q2**2 + q3**2

    return np.array([[R00, R01, R02], [R10, R11, R12], [R20, R21, R22]])


# --- 2. Set up the Scene ---
plotter = pv.Plotter(window_size=[1200, 800])

# Add static world axes
plotter.add_axes(
    line_width=3, xlabel="World X", ylabel="World Y", zlabel="World Z", color="black"
)

# --- 3. Create the Actors ---

# Create a cube mesh
cube = pv.Cube(center=(0, 0, 0), x_length=1.0, y_length=1.0, z_length=1.0)
cube_actor = plotter.add_mesh(
    cube,
    color="cornflowerblue",
    opacity=0.6,
    show_edges=True,
    edge_color="black",
    line_width=2,
)

# Create the body-fixed frame (three separate arrows)
arrow_x = pv.Arrow(
    direction=(1, 0, 0), tip_length=0.25, tip_radius=0.1, shaft_radius=0.03, scale=1.0
)
arrow_y = pv.Arrow(
    direction=(0, 1, 0), tip_length=0.25, tip_radius=0.1, shaft_radius=0.03, scale=1.0
)
arrow_z = pv.Arrow(
    direction=(0, 0, 1), tip_length=0.25, tip_radius=0.1, shaft_radius=0.03, scale=1.0
)

body_axes_actor_x = plotter.add_mesh(arrow_x, color="red", line_width=5)
body_axes_actor_y = plotter.add_mesh(arrow_y, color="green", line_width=5)
body_axes_actor_z = plotter.add_mesh(arrow_z, color="blue", line_width=5)


# --- 4. Define the Callback Function ---
def update_rotation(value):
    # This check prevents errors during initialization
    if not hasattr(plotter, "slider_widgets") or len(plotter.slider_widgets) < 4:
        return

    # Get the 4 slider values
    phi_deg = plotter.slider_widgets[0].GetRepresentation().GetValue()
    e_x = plotter.slider_widgets[1].GetRepresentation().GetValue()
    e_y = plotter.slider_widgets[2].GetRepresentation().GetValue()
    e_z = plotter.slider_widgets[3].GetRepresentation().GetValue()

    # --- Math Calculations ---

    # 1. Normalize the axis vector
    axis_vec = np.array([e_x, e_y, e_z])
    norm = np.linalg.norm(axis_vec)

    if norm < 1e-6:
        # Undefined axis (0,0,0) = no rotation
        e_normalized = np.array([1.0, 0.0, 0.0])  # Show a default axis
        q = np.array([1.0, 0.0, 0.0, 0.0])
        R = np.identity(3)
    else:
        e_normalized = axis_vec / norm

        # 2. Calculate Euler Parameters (q) from Axis-Angle
        phi_rad = np.radians(phi_deg)
        c = np.cos(phi_rad / 2.0)
        s = np.sin(phi_rad / 2.0)

        q0 = c
        q1 = e_normalized[0] * s
        q2 = e_normalized[1] * s
        q3 = e_normalized[2] * s
        q = np.array([q0, q1, q2, q3])

        # 3. Calculate DCM [R] from Euler Parameters
        R = quaternion_to_dcm(q)

    # 4. Create 4x4 transform matrix
    transform_matrix = np.identity(4)
    transform_matrix[:3, :3] = R

    # 5. Apply the transform to the rotating actors
    cube_actor.user_matrix = transform_matrix
    body_axes_actor_x.user_matrix = transform_matrix
    body_axes_actor_y.user_matrix = transform_matrix
    body_axes_actor_z.user_matrix = transform_matrix

    # 6. Format and update the text
    dcm_string = (
        "Direction Cosine Matrix [R]:\n"
        f"[{R[0,0]: 7.2f} {R[0,1]: 7.2f} {R[0,2]: 7.2f}]\n"
        f"[{R[1,0]: 7.2f} {R[1,1]: 7.2f} {R[1,2]: 7.2f}]\n"
        f"[{R[2,0]: 7.2f} {R[2,1]: 7.2f} {R[2,2]: 7.2f}]"
    )

    q_string = (
        "Euler Parameters (q):\n"
        f"q0 (scalar) = {q[0]: 7.3f}\n"
        f"q1 (vector) = {q[1]: 7.3f}\n"
        f"q2 (vector) = {q[2]: 7.3f}\n"
        f"q3 (vector) = {q[3]: 7.3f}"
    )

    # Combine the strings
    info_string = dcm_string + "\n\n" + q_string

    # 7. Create the (non-rotating) axis of rotation arrow
    axis_arrow_mesh = pv.Arrow(
        direction=e_normalized,
        scale=2.0,  # Make it long enough to see
        tip_length=0.2,
        tip_radius=0.08,
        shaft_radius=0.02,
    )

    # 8. Remove-and-replace logic for text and axis actors
    try:
        plotter.remove_actor(plotter._title_text_actor)
        plotter.remove_actor(plotter._info_text_actor)
        plotter.remove_actor(plotter._axis_actor)
    except AttributeError:
        pass  # First run

    plotter._title_text_actor = plotter.add_text(
        "Euler Parameters (from Axis-Angle)", font_size=15, position="upper_left"
    )

    # Add the combined text block
    plotter._info_text_actor = plotter.add_text(
        info_string, position="upper_right", font_size=12, color="black"
    )

    # Add the new axis arrow
    plotter._axis_actor = plotter.add_mesh(axis_arrow_mesh, color="purple", opacity=0.8)


# --- 5. Add the Sliders ---
# *** REPLACE event_type='always' with interaction_event='always' ***

plotter.add_slider_widget(
    callback=update_rotation,
    rng=[-180, 180],
    value=0.0,
    title="Angle (Φ)",
    pointa=(0.7, 0.4),
    pointb=(0.9, 0.4),
    interaction_event="always",
)
plotter.add_slider_widget(
    callback=update_rotation,
    rng=[-1, 1],
    value=1.0,
    title="Axis e_x",
    pointa=(0.7, 0.3),
    pointb=(0.9, 0.3),
    interaction_event="always",
)
plotter.add_slider_widget(
    callback=update_rotation,
    rng=[-1, 1],
    value=0.0,
    title="Axis e_y",
    pointa=(0.7, 0.2),
    pointb=(0.9, 0.2),
    interaction_event="always",
)
plotter.add_slider_widget(
    callback=update_rotation,
    rng=[-1, 1],
    value=0.0,
    title="Axis e_z",
    pointa=(0.7, 0.1),
    pointb=(0.9, 0.1),
    interaction_event="always",
)

# --- 6. Show the Plot ---
plotter.enable_parallel_projection()
plotter.enable_lightkit()

# Call update_rotation once manually to initialize the text
update_rotation(0.0)

plotter.show()
