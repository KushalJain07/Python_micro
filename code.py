import pandas as pd
from scipy import stats
from Bio import SeqIO
import numpy as np
import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

def convert(fasta_file):
    """
    Reads a fasta file, calculates sequence lengths, and returns a dataframe and stats.
    """

    try:
        records = SeqIO.parse(fasta_file, "fasta")
    except FileNotFoundError:
        messagebox.showerror("Error", "File is not available.")
        return None, None
    except Exception as e:
        messagebox.showerror("Error", f"Error reading fasta file: {e}")
        return None, None

    ids = []
    sequences = []
    lengths = []

    for record in records:
        ids.append(record.id)
        sequences.append(str(record.seq))
        lengths.append(len(record.seq))

    if not ids:
        messagebox.showwarning("Warning", "No sequences found in FASTA file.")
        return None, None

    df = pd.DataFrame({
        'Id': ids,
        'Sequence': sequences,
        'Length': lengths
    })

    # Calculate statistics
    mean_len = df['Length'].mean()
    median_len = df['Length'].median()
    mode_result = stats.mode(df['Length'], keepdims=True)

    stats_dict = {
        "mean": mean_len,
        "median": median_len,
        "mode": mode_result.mode[0],
        "mode_count": mode_result.count[0]
    }

    # Save CSV output
    output_filename = "project_output.csv"
    df.to_csv(output_filename, index=False)

    return df, stats_dict


# ---------------- GUI ---------------- #

def browse_file():
    filename = filedialog.askopenfilename(
        title="Select FASTA File",
        filetypes=[("FASTA Files", "*.fasta *.fa"), ("All Files", "*.*")]
    )
    if filename:
        file_path_var.set(filename)
        process_file()


def process_file():
    fasta = file_path_var.get().strip()

    if not fasta:
        messagebox.showwarning("Warning", "Please select a FASTA file first.")
        return

    df, stats_dict = convert(fasta)

    if df is not None:
        show_table(df)
        show_stats(stats_dict)

        messagebox.showinfo("Success", "Conversion complete.\nCSV saved as project_output.csv")


def show_table(df):
    for widget in table_frame.winfo_children():
        widget.destroy()

    tree = ttk.Treeview(table_frame)
    tree.pack(fill=tk.BOTH, expand=True)

    # Define columns
    tree["columns"] = list(df.columns)
    tree["show"] = "headings"

    for col in df.columns:
        tree.heading(col, text=col)
        tree.column(col, width=200)

    # Insert rows
    for _, row in df.iterrows():
        tree.insert("", tk.END, values=list(row))


def show_stats(stats):
    stat_str = (
        f"Mean Length           : {stats['mean']}\n"
        f"Median Length         : {stats['median']}\n"
        f"Mode Length           : {stats['mode']} "
        f"(appeared {stats['mode_count']} times)"
    )
    stats_label.config(text=stat_str)


# ---------------- MAIN WINDOW ---------------- #

root = tk.Tk()
root.title("FASTA File Analyzer")
root.geometry("900x600")

file_path_var = tk.StringVar()

# Top frame
top_frame = tk.Frame(root)
top_frame.pack(fill=tk.X, pady=10)

tk.Label(top_frame, text="Selected File: ").pack(side=tk.LEFT, padx=5)
tk.Entry(top_frame, textvariable=file_path_var, width=60).pack(side=tk.LEFT)
tk.Button(top_frame, text="Browse", command=browse_file).pack(side=tk.LEFT, padx=5)

# Table frame
table_frame = tk.Frame(root, bd=2, relief=tk.SUNKEN)
table_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

# Stats display
stats_label = tk.Label(root, font=("Courier", 10), justify=tk.LEFT)
stats_label.pack(pady=10)

root.mainloop()
