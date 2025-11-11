# Euler 321 Rotation Demo
# Vibe coded via Gemini

import pyvista as pv
import numpy as np

# --- 2. Set up the Scene (Moved to top) ---
plotter = pv.Plotter(window_size=[1200, 800])

# Add static world axes
plotter.add_axes(
    line_width=3, xlabel="World X", ylabel="World Y", zlabel="World Z", color="black"
)

# --- 3. Create the Actors (Moved up) ---

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


# --- 1. Define the Callback Function (Now can see actors) ---
def update_rotation(value):
    # Get the current values from the sliders
    # This check prevents the IndexError warnings
    if not hasattr(plotter, "slider_widgets") or len(plotter.slider_widgets) < 3:
        return  # Sliders are not ready yet

    yaw = plotter.slider_widgets[0].GetRepresentation().GetValue()
    pitch = plotter.slider_widgets[1].GetRepresentation().GetValue()
    roll = plotter.slider_widgets[2].GetRepresentation().GetValue()

    # Create a transformation object
    transform = pv.Transform()

    # Apply the Z-Y-X rotation sequence (Yaw, Pitch, Roll)
    transform.rotate_z(yaw)
    transform.rotate_y(pitch)
    transform.rotate_x(roll)

    # Apply the final 4x4 transformation matrix to the actors
    cube_actor.user_matrix = transform.matrix
    body_axes_actor_x.user_matrix = transform.matrix
    body_axes_actor_y.user_matrix = transform.matrix
    body_axes_actor_z.user_matrix = transform.matrix

    # Get the 3x3 rotation matrix (DCM)
    R = transform.matrix[:3, :3]

    # Format it into a clean string
    dcm_string = (
        "Direction Cosine Matrix (R):\n"
        f"[{R[0,0]: 7.2f} {R[0,1]: 7.2f} {R[0,2]: 7.2f}]\n"
        f"[{R[1,0]: 7.2f} {R[1,1]: 7.2f} {R[1,2]: 7.2f}]\n"
        f"[{R[2,0]: 7.2f} {R[2,1]: 7.2f} {R[2,2]: 7.2f}]"
    )

    # *** THIS IS THE NEW, CORRECTED LOGIC ***

    # 1. Try to remove the *previous* text actors.
    #    We use the private variable names (with '_')
    try:
        plotter.remove_actor(plotter._title_text_actor)
        plotter.remove_actor(plotter._dcm_text_actor)
    except AttributeError:
        # This happens on the first run, no actors to remove yet
        pass

    # 2. Add the new text actors and *store them*
    #    as private variables (with '_')
    plotter._title_text_actor = plotter.add_text(
        "Z-Y-X (Yaw-Pitch-Roll) Rotation", font_size=15, position="upper_left"
    )

    plotter._dcm_text_actor = plotter.add_text(
        dcm_string, position="upper_right", font_size=12, color="black"
    )


# --- 4. Add the Sliders (Now that callback is defined) ---

plotter.add_slider_widget(
    callback=update_rotation,
    rng=[-180, 180],
    value=0.0,
    title="Yaw (ψ)",
    pointa=(0.75, 0.3),
    pointb=(0.95, 0.3),
)
plotter.add_slider_widget(
    callback=update_rotation,
    rng=[-180, 180],
    value=0.0,
    title="Pitch (θ)",
    pointa=(0.75, 0.2),
    pointb=(0.95, 0.2),
)
plotter.add_slider_widget(
    callback=update_rotation,
    rng=[-180, 180],
    value=0.0,
    title="Roll (ϕ)",
    pointa=(0.75, 0.1),
    pointb=(0.95, 0.1),
)

# --- 6. Show the Plot (Final Step) ---
plotter.enable_parallel_projection()
plotter.enable_lightkit()

# Call update_rotation once manually to initialize the text
update_rotation(0.0)

plotter.show()
