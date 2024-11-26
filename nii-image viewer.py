import numpy as np
import tkinter as tk
from tkinter import filedialog, messagebox
import SimpleITK as sitk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class ImageViewer:
    def __init__(self, master):
        self.master = master
        self.master.title('nii Image Viewer')
        self.master.geometry("1200x600")
        self.top_frame = tk.Frame(master)
        self.left_frame = tk.Frame(master)
        self.right_frame = tk.Frame(master)
        self.bottom_frame = tk.Frame(master)
        self.top_frame.grid(row=0, column=0, columnspan=2, sticky="nsew")
        self.left_frame.grid(row=1, column=0, sticky="nsew")
        self.right_frame.grid(row=1, column=1, sticky="nsew")
        self.bottom_frame.grid(row=2, column=0, columnspan=2, sticky="nsew")
        self.master.grid_rowconfigure(1, weight=1)
        self.master.grid_columnconfigure(1, weight=1)
        self.load_button = tk.Button(self.top_frame, text="Load Image", command=self.load_image)
        self.load_button.pack(side=tk.TOP, pady=10)
        self.volume_label = tk.Label(self.top_frame, text="")
        self.volume_label.pack(side=tk.TOP, pady=10)
        self.fig, (self.ax1, self.ax2, self.ax3) = plt.subplots(1, 3, figsize=(15, 5))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.right_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.image = None
        self.slices = [0, 0, 0]
        self.image_info_label = tk.Label(self.top_frame, text="")
        self.image_info_label.pack(side=tk.TOP)
        self.create_left_panel()
        self.create_view_controls()
        self.cid_click = [self.fig.canvas.mpl_connect('button_press_event', self.onclick) for _ in range(3)]
        self.cid_move = [self.fig.canvas.mpl_connect('motion_notify_event', self.onmove) for _ in range(3)]
        self.cid_release = [self.fig.canvas.mpl_connect('button_release_event', self.onrelease) for _ in range(3)]
        self.click_coords = [None, None, None]
        self.measure_points = []
        self.measure_lines = []
        self.measure_polygons = []
        self.measure_labels = []
        self.measure_mode = None
        self.drawing = False
        self.last_x = None
        self.last_y = None
        self.current_polygon = []
        self.zoom_factor = 1.0
        self.cid_scroll = self.fig.canvas.mpl_connect('scroll_event', self.on_scroll)

    def create_left_panel(self):
        left_panel = self.left_frame

        self.window_center_label = tk.Label(left_panel, text="Window Center:")
        self.window_center_label.pack(pady=5)
        self.window_center_entry = tk.Entry(left_panel, width=10)
        self.window_center_entry.insert(0, "0")
        self.window_center_entry.pack(pady=5)

        self.window_width_label = tk.Label(left_panel, text="Window Width:")
        self.window_width_label.pack(pady=5)
        self.window_width_entry = tk.Entry(left_panel, width=10)
        self.window_width_entry.insert(0, "2000")
        self.window_width_entry.pack(pady=5)

        update_window_level_button = tk.Button(left_panel, text="Update Window Level", command=self.update_window_level)
        update_window_level_button.pack(pady=5)

        self.measure_length_button = tk.Button(left_panel, text="Measure Length", command=self.toggle_measure_length)
        self.measure_length_button.pack(pady=5)

        self.measure_area_button = tk.Button(left_panel, text="Measure Area", command=self.start_measure_area)
        self.measure_area_button.pack(pady=5)

        self.erase_button = tk.Button(left_panel, text="Erase", command=self.start_erase)
        self.erase_button.pack(pady=5)
    def create_view_controls(self):
        control_frame = self.bottom_frame
        labels = ['Axial', 'Coronal', 'Sagittal']
        self.slider_axial = None
        self.slider_coronal = None
        self.slider_sagittal = None

        for i, label in enumerate(labels):
            view_frame = tk.Frame(control_frame)
            view_frame.pack(side=tk.LEFT, fill=tk.X, expand=True)
            slider = tk.Scale(view_frame, from_=0, to=100, orient=tk.HORIZONTAL, label=f'{label} Slice',
                              command=lambda value, idx=i: self.update_slice(value, idx))
            slider.pack(fill=tk.X, expand=True)
            label = tk.Label(view_frame, text=label)
            label.pack(fill=tk.X)

            if i == 0:
                self.slider_axial = slider
            elif i == 1:
                self.slider_coronal = slider
            elif i == 2:
                self.slider_sagittal = slider

    def load_image(self):
        filename = filedialog.askopenfilename(title="Select Image",
                                              filetypes=(("NIfTI files", "*.nii"), ("all files", "*.*")))
        if filename:
            self.image = sitk.ReadImage(filename)
            self.image_array = sitk.GetArrayFromImage(self.image)
            self.slices = [self.image_array.shape[0] // 2, self.image_array.shape[1] // 2,
                           self.image_array.shape[2] // 2]

            self.slider_axial.config(to=self.image_array.shape[0] - 1)
            self.slider_coronal.config(to=self.image_array.shape[1] - 1)
            self.slider_sagittal.config(to=self.image_array.shape[2] - 1)

            self.slider_axial.set(self.slices[0])
            self.slider_coronal.set(self.slices[1])
            self.slider_sagittal.set(self.slices[2])

            self.show_image_info()
            self.display_slices()

    def show_image_info(self):
        size = self.image.GetSize()
        spacing = self.image.GetSpacing()
        origin = self.image.GetOrigin()
        info_text = f"Size: {size}\nSpacing: {spacing}\nOrigin: {origin}"
        self.image_info_label.config(text=info_text)

    def update_slice(self, value, index):
        if self.image is not None:
            self.slices[index] = int(value)
            self.display_slices()

            if index == 0:
                self.slider_coronal.set(self.slices[1])
                self.slider_sagittal.set(self.slices[2])
            elif index == 1:
                self.slider_axial.set(self.slices[0])
                self.slider_sagittal.set(self.slices[2])
            elif index == 2:
                self.slider_axial.set(self.slices[0])
                self.slider_coronal.set(self.slices[1])

    def display_slices(self):
        self.ax1.cla()
        self.ax2.cla()
        self.ax3.cla()
        windowed_image = self.apply_window_level(self.image_array)
        self.ax1.imshow(windowed_image[self.slices[0], :, :], cmap='gray')
        self.ax2.imshow(windowed_image[:, self.slices[1], :], cmap='gray')
        self.ax3.imshow(windowed_image[:, :, self.slices[2]], cmap='gray')

        if self.click_coords[0] is not None:
            self.ax1.scatter(self.click_coords[2], self.click_coords[1], color='brown')
            self.ax2.scatter(self.click_coords[0], self.click_coords[2], color='brown')
            self.ax3.scatter(self.click_coords[0], self.click_coords[1], color='brown')

        for line in self.measure_lines:
            line[0].plot(line[1], line[2], color='red')
        for polygon in self.measure_polygons:
            polygon[0].fill(polygon[1], polygon[2], alpha=0.3, color='lightgreen')
        for label in self.measure_labels:
            label[0].text(label[1], label[2], label[3], fontsize=12, color='green')
        self.canvas.draw()
    def apply_window_level(self, image_array):
        try:
            self.window_center = float(self.window_center_entry.get())
            self.window_width = float(self.window_width_entry.get())
        except ValueError:
            pass

        min_value = self.window_center - self.window_width / 2
        max_value = self.window_center + self.window_width / 2
        windowed_image = np.clip(image_array, min_value, max_value)
        windowed_image = ((windowed_image - min_value) / (max_value - min_value)) * 255
        return windowed_image.astype(np.uint8)

    def update_window_level(self):
        self.display_slices()

    def onclick(self, event):
        if event.inaxes in [self.ax1, self.ax2, self.ax3]:
            ax_index = [self.ax1, self.ax2, self.ax3].index(event.inaxes)

            if event.xdata is not None and event.ydata is not None:
                x, y = int(event.xdata), int(event.ydata)
                i, j, k = -1, -1, -1
                if ax_index == 0:
                    i = self.slices[0]
                    j, k = y, x
                elif ax_index == 1:
                    i, k = y, self.slices[1]
                    j = x
                elif ax_index == 2:
                    i, j = self.slices[2], y
                    k = x

                self.slices = [i, j, k]
                self.update_slice(i, ax_index)
                self.click_coords = [i, j, k]

                point = self.image.TransformIndexToPhysicalPoint((k, j, i))
                self.image_info_label.config(
                    text=f"Size: {self.image.GetSize()}\nSpacing: {self.image.GetSpacing()}\nOrigin: {self.image.GetOrigin()}\nClicked at: (i, j, k) = ({i}, {j}, {k}), (x, y, z) = ({point[0]:.2f}, {point[1]:.2f}, {point[2]:.2f})")

                if self.measure_mode == 'area':
                    self.drawing = True
                    self.last_x, self.last_y = x, y
                    self.current_polygon.append((x, y))
                elif self.measure_mode == 'length':
                    self.measure_points.append((ax_index, x, y))
                    if len(self.measure_points) == 2:
                        self.measure_length()

    def onmove(self, event):
        if event.inaxes in [self.ax1, self.ax2, self.ax3] and self.drawing:
            ax_index = [self.ax1, self.ax2, self.ax3].index(event.inaxes)
            if event.xdata is not None and event.ydata is not None:
                x, y = int(event.xdata), int(event.ydata)
                if self.last_x is not None and self.last_y is not None:
                    ax = [self.ax1, self.ax2, self.ax3][ax_index]
                    ax.plot([self.last_x, x], [self.last_y, y], color='lightgreen')
                    self.canvas.draw()
                self.last_x, self.last_y = x, y
                self.current_polygon.append((x, y))

    def onrelease(self, event):
        if self.drawing:
            self.drawing = False
            self.measure_polygons.append(
                (self.ax1, [p[0] for p in self.current_polygon], [p[1] for p in self.current_polygon]))
            area = self.calculate_polygon_area(self.current_polygon)
            self.measure_labels.append(
                (self.ax1, self.current_polygon[-1][0], self.current_polygon[-1][1], f"Area: {area:.2f}"))
            self.current_polygon = []
            self.display_slices()

    def start_measure_area(self):
        self.measure_points.clear()
        self.measure_lines.clear()
        self.measure_labels.clear()
        self.measure_length_button.config(state=tk.NORMAL)
        self.measure_area_button.config(state=tk.DISABLED)
        self.measure_mode = 'area'

    def start_erase(self):
        self.measure_points.clear()
        self.measure_lines.clear()
        self.measure_polygons.clear()
        self.measure_labels.clear()
        self.measure_mode = None
        self.measure_length_button.config(state=tk.NORMAL)
        self.measure_area_button.config(state=tk.NORMAL)
        self.display_slices()

    def toggle_measure_length(self):
        if self.measure_mode != 'length':
            self.measure_points.clear()
            self.measure_lines.clear()
            self.measure_labels.clear()
            self.measure_length_button.config(state=tk.DISABLED)
            self.measure_area_button.config(state=tk.NORMAL)
            self.measure_mode = 'length'
        else:
            self.measure_mode = None
            self.measure_length_button.config(state=tk.NORMAL)

    def measure_length(self):
        self.measure_lines.clear()
        self.measure_labels.clear()

        point1 = self.measure_points[0][1:]
        point2 = self.measure_points[1][1:]
        distance = np.sqrt((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2)

        ax_index = self.measure_points[0][0]
        ax = [self.ax1, self.ax2, self.ax3][ax_index]
        line = ax.plot([point1[0], point2[0]], [point1[1], point2[1]], color='red')
        self.measure_lines.append((line[0], [point1[0], point2[0]], [point1[1], point2[1]]))

        label = ax.text((point1[0] + point2[0]) / 2, (point1[1] + point2[1]) / 2, f'{distance:.2f}', fontsize=12,
                        color='green')
        self.measure_labels.append((label, (point1[0] + point2[0]) / 2, (point1[1] + point2[1]) / 2, f'{distance:.2f}'))

        self.canvas.draw()
        self.measure_mode = None
        self.measure_length_button.config(state=tk.NORMAL)

    def calculate_polygon_area(self, polygon):
        x = [p[0] for p in polygon]
        y = [p[1] for p in polygon]
        return 0.5 * abs(sum(x[i] * y[i - 1] - x[i - 1] * y[i] for i in range(len(polygon))))
    def on_scroll(self, event):
        if event.inaxes in [self.ax1, self.ax2, self.ax3]:
            zoom_in = event.button == 'up'
            factor = 1.2 if zoom_in else 1/1.2
            self.zoom_factor *= factor
            ax = event.inaxes
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            new_xlim = [(x - event.xdata) * factor + event.xdata for x in xlim]
            new_ylim = [(y - event.ydata) * factor + event.ydata for y in ylim]
            ax.set_xlim(new_xlim)
            ax.set_ylim(new_ylim)
            self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = ImageViewer(root)
    root.mainloop()