import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog
from tkinter import ttk

def complement_sequence(sequence):
    dict1 = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    complement = ""
    for i in sequence:
        reverse = dict1.get(i)
        if reverse is not None:
            complement += reverse
        else:
            complement += 'N'
    return complement

def transcribe_to_rna(dna_sequence):
    return dna_sequence.replace('T', 'U')

def reverse_transcribe_to_dna(rna_sequence):
    return rna_sequence.replace('U', 'T')

dna_codon_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'STOP', 'TAG': 'STOP',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'STOP', 'TGG': 'W'
}

rna_codon_table = {
    'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
    'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
    'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
    'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
    'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
    'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
    'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
    'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
    'UAC': 'Y', 'UAU': 'Y', 'UAA': 'STOP', 'UAG': 'STOP',
    'UGC': 'C', 'UGU': 'C', 'UGA': 'STOP', 'UGG': 'W'
}

def translate_sequence(sequence, start_index, is_rna=False):
    codon_table = rna_codon_table if is_rna else dna_codon_table
    if 'T' in sequence and 'U' in sequence:
        raise ValueError("Invalid sequence. Both 'T' and 'U' found in the same sequence.")
    
    if is_rna:
        sequence = reverse_transcribe_to_dna(sequence)
    sequence_length = len(sequence)
    protein_sequence = ""
    start_codon_found = False

    for i in range(start_index, sequence_length, 3):
        codon = sequence[i:i + 3]
        if codon == 'ATG':
            start_codon_found = True
        if start_codon_found:
            if codon in ['TAA', 'TAG', 'TGA']:
                break
            amino_acid = codon_table.get(codon, 'X')
            protein_sequence += amino_acid

    return protein_sequence if start_codon_found else None

def translate_reverse_reading_frames(sequence, is_rna=False):
    complement_seq = complement_sequence(sequence[::-1])

    for frame in range(3):
        frame_text = f"Reverse Complement Reading Frame {frame + 1}:"
        protein_sequence = translate_sequence(complement_seq, frame, is_rna)
        if protein_sequence:
            frame_text += f"\nProtein sequence: {protein_sequence}"
        else:
            frame_text += "\nNo START codon found"
        results_text.insert(tk.END, frame_text + "\n\n")

def handle_translation():
    results_text.delete('1.0', tk.END)
    input_type = input_var.get().upper()

    if input_type not in ('DNA', 'RNA'):
        messagebox.showerror("Error", "Invalid input. Please enter 'DNA' or 'RNA'.")
        return
    input_sequence = sequence_entry.get().upper()

    if not all(base in 'ATGCU' for base in input_sequence):
        messagebox.showerror("Error", "Invalid sequence. Please enter a valid RNA or DNA sequence.")
        return

    try:
        is_rna = input_type == 'RNA'
        for frame in range(3):
            frame_text = f"Reading Frame {frame + 1}:"
            protein_sequence = translate_sequence(input_sequence, frame, is_rna)
            if protein_sequence:
                frame_text += f"\nProtein sequence: {protein_sequence}"
            else:
                frame_text += "\nNo START codon found"
            results_text.insert(tk.END, frame_text + "\n\n")
        translate_reverse_reading_frames(input_sequence, is_rna)
    except ValueError as e:
        messagebox.showerror("Error", str(e))

def open_file():
    file_path = filedialog.askopenfilename(filetypes=[("Text files", "*.txt")])
    if file_path:
        try:
            with open(file_path, 'r') as file:
                sequence = file.read().strip()
                sequence_entry.delete(0, tk.END)
                sequence_entry.insert(0, sequence)
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred while opening the file: {e}")

def reset_fields():
    input_var.set('')
    sequence_entry.delete(0, tk.END)
    results_text.delete('1.0', tk.END)

# Tooltip function
def create_tooltip(widget, text):
    tool_tip = ttk.Label(root, text=text)
    tool_tip.place(in_=widget, relx=1.0, rely=0, x=-10, y=1, anchor="ne")
    widget.bind("<Enter>", lambda event: tool_tip.place_configure(relx=1.0, rely=0, x=-10, y=1, anchor="ne"))
    widget.bind("<Leave>", lambda event: tool_tip.place_forget())

# GUI Setup
root = tk.Tk()
root.title("DNA/RNA Translation")

root.title("DNA/RNA Translation")
root.configure(bg='#d3e0ea')  # Set background color

header_label = tk.Label(root, text="DNA/RNA Translation", font=("Arial", 18, "bold"), bg='#f0f0f0', fg='#333')
header_label.pack(pady=(20, 10))

input_frame = tk.Frame(root, bg='#d3e0ea')  # Set background color
input_frame.pack(pady=(10, 20))

input_label = tk.Label(input_frame, text="Enter 'DNA' or 'RNA':", bg='#f0f0f0', font=("Arial", 12), fg='#333')
input_label.grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)

input_var = tk.StringVar()
input_entry = tk.Entry(input_frame, textvariable=input_var, font=("Arial", 12))
input_entry.grid(row=0, column=1, padx=5, pady=5)
create_tooltip(input_entry, "Enter 'DNA' or 'RNA'")

sequence_label = tk.Label(input_frame, text="Enter sequence:", bg='#f0f0f0', font=("Arial", 12), fg='#333')
sequence_label.grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)

sequence_entry = tk.Entry(input_frame, font=("Arial", 12))
sequence_entry.grid(row=1, column=1, padx=5, pady=5)
create_tooltip(sequence_entry, "Enter DNA or RNA sequence")

file_button = tk.Button(input_frame, text="Open File", command=open_file, bg='#4caf50', fg='white', font=("Arial", 12, "bold"))
file_button.grid(row=2, column=0, columnspan=2, padx=5, pady=15)

translate_button = tk.Button(input_frame, text="Translate", command=handle_translation, bg='#4caf50', fg='white', font=("Arial", 12, "bold"))
translate_button.grid(row=3, column=0, columnspan=2, padx=5, pady=15)

reset_button = tk.Button(input_frame, text="Reset", command=reset_fields, bg='#f44336', fg='white', font=("Arial", 12, "bold"))
reset_button.grid(row=4, column=0, columnspan=2, padx=5, pady=15)

# Add scrollbar
scrollbar = tk.Scrollbar(root)
scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

results_text = tk.Text(root, height=15, width=60, font=("Arial", 12), yscrollcommand=scrollbar.set)
results_text.pack(pady=(0, 20))
scrollbar.config(command=results_text.yview)

root.mainloop()


