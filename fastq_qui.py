import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from tkinterdnd2 import DND_FILES, TkinterDnD
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import threading
import os
from fastaq import FastqReader

class FastqAnalyzerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("FASTQ File Analyzer")
        self.root.geometry("1000x700")
        
        self.current_file = None
        self.analyzer = None
        
        self.setup_gui()
        
    def setup_gui(self):
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        
        # File selection section
        file_frame = ttk.LabelFrame(main_frame, text="File Selection", padding="10")
        file_frame.grid(row=0, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))
        file_frame.columnconfigure(1, weight=1)
        
        ttk.Button(file_frame, text="Browse FASTQ File", 
                  command=self.browse_file).grid(row=0, column=0, padx=(0, 10))
        
        self.file_label = ttk.Label(file_frame, text="No file selected")
        self.file_label.grid(row=0, column=1, sticky=(tk.W, tk.E))
        
        # Drop area for drag and drop
        drop_frame = ttk.Frame(file_frame, relief="sunken", height=50)
        drop_frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(10, 0))
        drop_frame.columnconfigure(0, weight=1)
        
        drop_label = ttk.Label(drop_frame, text="Or drag and drop FASTQ file here", 
                              foreground="gray")
        drop_label.grid(row=0, column=0, pady=15)
        
        # Register for drag and drop
        drop_frame.drop_target_register(DND_FILES)
        drop_frame.dnd_bind('<<Drop>>', self.on_drop)
        
        # Analysis controls
        control_frame = ttk.LabelFrame(main_frame, text="Analysis Controls", padding="10")
        control_frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))
        
        ttk.Button(control_frame, text="Show Statistics", 
                  command=self.show_statistics).grid(row=0, column=0, padx=(0, 10))
        
        ttk.Button(control_frame, text="Quality Analysis", 
                  command=self.quality_analysis).grid(row=0, column=1, padx=(0, 10))
        
        ttk.Button(control_frame, text="Nucleotide Content", 
                  command=self.nucleotide_analysis).grid(row=0, column=2, padx=(0, 10))
        
        ttk.Button(control_frame, text="Full Analysis", 
                  command=self.full_analysis).grid(row=0, column=3)
        
        # Progress bar
        self.progress = ttk.Progressbar(control_frame, mode='indeterminate')
        self.progress.grid(row=1, column=0, columnspan=4, sticky=(tk.W, tk.E), pady=(10, 0))
        
        # Results section
        results_frame = ttk.LabelFrame(main_frame, text="Analysis Results", padding="10")
        results_frame.grid(row=2, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(0, 10))
        results_frame.columnconfigure(0, weight=1)
        results_frame.rowconfigure(0, weight=1)
        
        # Create notebook for tabs
        self.notebook = ttk.Notebook(results_frame)
        self.notebook.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Statistics tab
        self.stats_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.stats_frame, text="Statistics")
        
        self.stats_text = tk.Text(self.stats_frame, wrap=tk.WORD, width=80, height=10)
        stats_scrollbar = ttk.Scrollbar(self.stats_frame, orient="vertical", command=self.stats_text.yview)
        self.stats_text.configure(yscrollcommand=stats_scrollbar.set)
        
        self.stats_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        stats_scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
        self.stats_frame.columnconfigure(0, weight=1)
        self.stats_frame.rowconfigure(0, weight=1)
        
        # Plots tab
        self.plots_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.plots_frame, text="Plots")
        
        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, relief="sunken")
        status_bar.grid(row=3, column=0, columnspan=2, sticky=(tk.W, tk.E))
        
        # Configure main frame grid weights
        main_frame.rowconfigure(2, weight=1)
        
    def browse_file(self):
        """Open file dialog to select FASTQ file"""
        filename = filedialog.askopenfilename(
            title="Select FASTQ File",
            filetypes=[
                ("FASTQ files", "*.fastq *.fq"),
                ("Compressed FASTQ", "*.fastq.gz *.fq.gz"),
                ("All files", "*.*")
            ]
        )
        
        if filename:
            self.load_file(filename)
    
    def on_drop(self, event):
        """Handle drag and drop file"""
        # Extract filename from drop event
        filename = event.data.strip('{}')  # Remove curly braces on Windows
        if os.path.isfile(filename):
            self.load_file(filename)
    
    def load_file(self, filename):
        """Load and validate FASTQ file"""
        self.current_file = filename
        self.file_label.config(text=os.path.basename(filename))
        self.status_var.set(f"Loaded: {os.path.basename(filename)}")
        
        # Clear previous results
        self.stats_text.delete(1.0, tk.END)
        self.clear_plots()
    
    def clear_plots(self):
        """Clear all plots from plots frame"""
        for widget in self.plots_frame.winfo_children():
            widget.destroy()
    
    def show_statistics(self):
        """Display basic statistics"""
        if not self.current_file:
            messagebox.showerror("Error", "Please select a FASTQ file first")
            return
        
        self.start_processing()
        
        def analyze():
            try:
                self.analyzer = FastqReader(self.current_file)
                count = self.analyzer.get_sequence_count()
                avg_len = self.analyzer.get_average_length()
                total_bp = count * avg_len
                
                # Update GUI in main thread
                self.root.after(0, lambda: self.display_statistics(count, avg_len, total_bp))
                
            except Exception as e:
                self.root.after(0, lambda: messagebox.showerror("Error", f"Analysis failed: {str(e)}"))
            finally:
                self.root.after(0, self.stop_processing)
        
        # Run analysis in separate thread
        thread = threading.Thread(target=analyze)
        thread.daemon = True
        thread.start()
    
    def display_statistics(self, count, avg_len, total_bp):
        """Display statistics in the text widget"""
        stats_text = f"""FASTQ FILE STATISTICS
File: {os.path.basename(self.current_file)}

Basic Statistics:
• Sequence count: {count:,}
• Average length: {avg_len:.2f} bp
• Total data volume: {total_bp:,.0f} bp

File Information:
• File size: {os.path.getsize(self.current_file) / 1024 / 1024:.2f} MB
• Estimated memory usage: Optimized (generators)

Analysis completed successfully."""
        
        self.stats_text.delete(1.0, tk.END)
        self.stats_text.insert(1.0, stats_text)
        self.notebook.select(0)  # Switch to statistics tab
    
    def quality_analysis(self):
        """Perform quality analysis and display plots"""
        if not self.current_file:
            messagebox.showerror("Error", "Please select a FASTQ file first")
            return
        
        self.start_processing()
        
        def analyze():
            try:
                if not self.analyzer:
                    self.analyzer = FastqReader(self.current_file)
                
                # Create plots
                self.analyzer.plot_per_base_quality("temp_quality.png")
                self.analyzer.plot_sequence_length_distribution("temp_length.png")
                
                # Update GUI in main thread
                self.root.after(0, self.display_quality_plots)
                
            except Exception as e:
                self.root.after(0, lambda: messagebox.showerror("Error", f"Quality analysis failed: {str(e)}"))
            finally:
                self.root.after(0, self.stop_processing)
        
        thread = threading.Thread(target=analyze)
        thread.daemon = True
        thread.start()
    
    def display_quality_plots(self):
        """Display quality plots in the GUI"""
        self.clear_plots()
        
        # Create frame for plots
        plots_container = ttk.Frame(self.plots_frame)
        plots_container.pack(fill=tk.BOTH, expand=True)
        
        try:
            # Display quality plot
            fig1 = plt.figure(figsize=(8, 4))
            img1 = plt.imread("temp_quality.png")
            plt.imshow(img1)
            plt.axis('off')
            plt.title('Per Base Sequence Quality', pad=20)
            
            canvas1 = FigureCanvasTkAgg(fig1, plots_container)
            canvas1.draw()
            canvas1.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)
            
            # Display length distribution plot
            fig2 = plt.figure(figsize=(8, 4))
            img2 = plt.imread("temp_length.png")
            plt.imshow(img2)
            plt.axis('off')
            plt.title('Sequence Length Distribution', pad=20)
            
            canvas2 = FigureCanvasTkAgg(fig2, plots_container)
            canvas2.draw()
            canvas2.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
            
            # Switch to plots tab
            self.notebook.select(1)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to display plots: {str(e)}")
    
    def nucleotide_analysis(self):
        """Perform nucleotide content analysis"""
        if not self.current_file:
            messagebox.showerror("Error", "Please select a FASTQ file first")
            return
        
        self.start_processing()
        
        def analyze():
            try:
                if not self.analyzer:
                    self.analyzer = FastqReader(self.current_file)
                
                self.analyzer.plot_per_base_content("temp_content.png")
                self.root.after(0, self.display_nucleotide_plot)
                
            except Exception as e:
                self.root.after(0, lambda: messagebox.showerror("Error", f"Content analysis failed: {str(e)}"))
            finally:
                self.root.after(0, self.stop_processing)
        
        thread = threading.Thread(target=analyze)
        thread.daemon = True
        thread.start()
    
    def display_nucleotide_plot(self):
        """Display nucleotide content plot"""
        self.clear_plots()
        
        try:
            fig = plt.figure(figsize=(10, 6))
            img = plt.imread("temp_content.png")
            plt.imshow(img)
            plt.axis('off')
            plt.title('Nucleotide Content Analysis', pad=20)
            
            canvas = FigureCanvasTkAgg(fig, self.plots_frame)
            canvas.draw()
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
            
            self.notebook.select(1)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to display nucleotide plot: {str(e)}")
    
    def full_analysis(self):
        """Perform complete analysis"""
        if not self.current_file:
            messagebox.showerror("Error", "Please select a FASTQ file first")
            return
        
        self.start_processing()
        
        def analyze():
            try:
                self.analyzer = FastqReader(self.current_file)
                
                # Get statistics
                count = self.analyzer.get_sequence_count()
                avg_len = self.analyzer.get_average_length()
                total_bp = count * avg_len
                
                # Create all plots
                self.analyzer.generate_all_plots()
                
                # Update GUI
                self.root.after(0, lambda: self.display_full_analysis(count, avg_len, total_bp))
                
            except Exception as e:
                self.root.after(0, lambda: messagebox.showerror("Error", f"Full analysis failed: {str(e)}"))
            finally:
                self.root.after(0, self.stop_processing)
        
        thread = threading.Thread(target=analyze)
        thread.daemon = True
        thread.start()
    
    def display_full_analysis(self, count, avg_len, total_bp):
        """Display results of full analysis"""
        # Show statistics
        self.display_statistics(count, avg_len, total_bp)
        
        # Show all plots
        self.display_quality_plots()
        
        messagebox.showinfo("Analysis Complete", "Full analysis completed successfully!\n\nCheck the Statistics and Plots tabs for results.")
    
    def start_processing(self):
        """Start progress indicator"""
        self.progress.start()
        self.status_var.set("Processing...")
    
    def stop_processing(self):
        """Stop progress indicator"""
        self.progress.stop()
        self.status_var.set("Ready")

def main():
    # Create and run the GUI application
    root = TkinterDnD.Tk()
    app = FastqAnalyzerGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
