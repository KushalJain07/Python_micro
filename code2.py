import traceback
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

# numeric & plotting libs
import numpy as np
import pandas as pd
from Bio import SeqIO

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


# ---------------- CORE ----------------

def convert(fasta_file):
    p = Path(fasta_file)
    

    with p.open("r") as fh:
        records = list(SeqIO.parse(fh, "fasta"))

    

    ids = [r.id for r in records]
    seqs = [str(r.seq) for r in records]
    lengths = [len(r.seq) for r in records]

    df = pd.DataFrame({"Id": ids, "Sequence": seqs, "Length": lengths})

    mean_len = float(df["Length"].mean())
    median_len = float(df["Length"].median())

    vc = df["Length"].value_counts().sort_values(ascending=False)
    mode_len = int(vc.index[0])
    mode_count = int(vc.iloc[0])

    stats = {
        "mean": mean_len,
        "median": median_len,
        "mode": mode_len,
        "mode_count": mode_count,
        "vc": vc
    }

    try:
        df.to_csv("project_output.csv", index=False)
    except:
        pass

    return df, stats


def make_bar_plot(stats):
    """Create wide & readable bar chart."""
    labels = ["Mean", "Median", "Mode"]
    values = [stats["mean"], stats["median"], stats["mode"]]

    fig = Figure(figsize=(6, 4), dpi=80)  # <- WIDER plot
    ax = fig.add_subplot(1, 1, 1)

    x = np.arange(len(labels))
    bars = ax.bar(x, values, edgecolor="black", linewidth=1)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_ylabel("Length", fontsize=12)
    ax.set_title("Mean | Median | Mode", fontsize=14)

    for rect, val in zip(bars, values):
        ax.annotate(f"{val:.2f}",
                     xy=(rect.get_x() + rect.get_width() / 2, rect.get_height()),
                     xytext=(0, 6), textcoords="offset points",
                     ha="center", fontsize=11, fontweight="bold")

    fig.tight_layout()
    return fig


# ---------------- GUI APP ----------------

class App:
    def __init__(self, root):
        self.root = root
        root.title("FASTA Analyzer — Table & Summary Graph")
        root.geometry("1350x750")

        self.current_fig = None

        # --- TOP BAR ---
        top = tk.Frame(root)
        top.pack(fill=tk.X, pady=10)

        tk.Label(top, text="FASTA File: ").pack(side=tk.LEFT, padx=5)
        self.file_var = tk.StringVar()
        tk.Entry(top, textvariable=self.file_var, width=60).pack(side=tk.LEFT)

        tk.Button(top, text="Browse", command=self.browse).pack(side=tk.LEFT, padx=10)
        tk.Button(top, text="Download Plot", command=self.download_plot).pack(side=tk.LEFT, padx=5)

        # --- SPLIT AREA ---
        main_split = tk.PanedWindow(root, orient=tk.HORIZONTAL)
        main_split.pack(fill=tk.BOTH, expand=True)

        # ------------------------------------------------------------------
        # --- MODIFICATIONS START HERE ---
        # ------------------------------------------------------------------

        # Wrapper Frame 1 (This goes into the PanedWindow)
        wrapper_left = tk.Frame(main_split)
        main_split.add(wrapper_left, stretch="always")

        # LEFT PANEL (TABLE) - This is the actual frame, placed inside the wrapper
        # The padx/pady here creates the 10-pixel margin/padding inside the wrapper.
        self.left = tk.Frame(wrapper_left, bd=2, relief=tk.SUNKEN)
        self.left.grid(row=0, column=0, sticky="NSEW", padx=10, pady=10)
        
        # Configure wrapper grid weight so self.left fills the space
        wrapper_left.grid_rowconfigure(0, weight=1)
        wrapper_left.grid_columnconfigure(0, weight=1)


        # Wrapper Frame 2
        wrapper_right = tk.Frame(main_split)
        main_split.add(wrapper_right, stretch="always")

        # RIGHT PANEL (GRAPH) - This is the actual frame, placed inside the wrapper
        # The padx/pady here creates the 10-pixel margin/padding inside the wrapper.
        self.right = tk.Frame(wrapper_right, bd=2, relief=tk.SUNKEN)
        self.right.grid(row=0, column=0, sticky="NSEW", padx=10, pady=10)

        # Configure wrapper grid weight so self.right fills the space
        wrapper_right.grid_rowconfigure(0, weight=1)
        wrapper_right.grid_columnconfigure(0, weight=1)
        
        # ------------------------------------------------------------------
        # --- MODIFICATIONS END HERE ---
        # ------------------------------------------------------------------


        # Set 50/50 division (Use a slight offset to account for margins)
        # Note: The PanedWindow now holds the wrappers, which internally manage the margin.
        root.after(100, lambda: main_split.sash_place(0, 650, 0))

        # --- TABLE SETUP ---
        # The table now needs to be packed into the self.left frame
        self.table = ttk.Treeview(self.left)
        self.table.pack(fill=tk.BOTH, expand=True)

        # Stats text under table
        # NOTE: You had two stats_label assignments here, fixed to use two labels for potential separate content
        self.stats_label_top = tk.Label(self.left, font=("Courier", 11), justify=tk.LEFT)
        self.stats_label_top.pack(pady=10)
        self.stats_label_bottom = tk.Label(self.left, font=("Courier", 11), justify=tk.RIGHT)
        self.stats_label_bottom.pack(pady=10)


    # ---------------- CALLBACKS ----------------

    def browse(self):
        path = filedialog.askopenfilename(
            title="Select FASTA File",
            filetypes=[("FASTA", "*.fasta *.fa *.fna"), ("All Files", "*.*")]
        )
        if path:
            self.file_var.set(path)
            self.process(path)

    def process(self, file):
        try:
            df, stats = convert(file)
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Error", str(e))
            return

        self.load_table(df)
        self.show_stats(stats)
        self.draw_plot(stats)

    def load_table(self, df):
        # Clear previous
        self.table.delete(*self.table.get_children())
        self.table["columns"] = list(df.columns)
        self.table["show"] = "headings"

        for col in df.columns:
            self.table.heading(col, text=col)
            self.table.column(col, width=180)

        for _, row in df.iterrows():
            self.table.insert("", tk.END, values=list(row))

    def show_stats(self, stats):
        txt = (
            f"Mean   : {stats['mean']:.2f}\n"
            f"Median : {stats['median']:.2f}\n"
            f"Mode   : {stats['mode']} (count={stats['mode_count']})"
        )
        # Using the top stats label for the main statistics
        self.stats_label_top.config(text=txt)

    def draw_plot(self, stats):
        for w in self.right.winfo_children():
            w.destroy()

        fig = make_bar_plot(stats)
        self.current_fig = fig

        canvas = FigureCanvasTkAgg(fig, master=self.right)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def download_plot(self):
        if not self.current_fig:
            messagebox.showwarning("No Plot", "No plot to save.")
            return

        name = filedialog.asksaveasfilename(defaultextension=".png",
                                             filetypes=[("PNG Image", "*.png")])
        if not name:
            return

        try:
            self.current_fig.savefig(name, dpi=300)
            messagebox.showinfo("Saved", "Plot saved successfully.")
        except Exception as e:
            messagebox.showerror("Error", str(e))


# ---------------- MAIN ----------------

def main():
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()