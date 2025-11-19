"""
FASTA Single-summary GUI (table left + wide graph right)

Creates ONE bar chart (Mean | Median | Mode)
and gives proper 50/50 layout without compacting the graph.
"""

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
    if not p.exists():
        raise FileNotFoundError(f"No file: {p}")

    with p.open("r") as fh:
        records = list(SeqIO.parse(fh, "fasta"))

    if not records:
        raise ValueError("No sequences found in FASTA")

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

    fig = Figure(figsize=(8, 5), dpi=120)   # <- WIDER plot
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

        # LEFT PANEL (TABLE)
        self.left = tk.Frame(main_split, bd=2, relief=tk.SUNKEN)
        main_split.add(self.left, stretch="always")

        # RIGHT PANEL (GRAPH)
        self.right = tk.Frame(main_split, bd=2, relief=tk.SUNKEN)
        main_split.add(self.right, stretch="always")

        # Set 50/50 division
        root.after(100, lambda: main_split.sash_place(0, 650, 0))

        # --- TABLE SETUP ---
        self.table = ttk.Treeview(self.left)
        self.table.pack(fill=tk.BOTH, expand=True)

        # Stats text under table
        self.stats_label = tk.Label(self.left, font=("Courier", 11), justify=tk.LEFT)
        self.stats_label.pack(pady=10)


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
        self.stats_label.config(text=txt)

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
