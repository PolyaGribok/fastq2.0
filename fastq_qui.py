import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from collections import defaultdict
import threading
import os

# Фисташковая цветовая схема
PISTACHIO_THEME = {
    'primary': '#93C572',      # Основной фисташковый
    'primary_light': '#B8D8A3', # Светлый фисташковый
    'primary_dark': '#7AA65A',  # Тёмный фисташковый
    'secondary': '#F5F5DC',     # Бежевый фон
    'accent': '#8FBC8F',        # Акцентный цвет
    'text_dark': '#2F4F2F',     # Тёмно-зелёный текст
    'text_light': '#FFFFFF',    # Белый текст
    'success': '#90EE90',       # Зелёный успех
    'warning': '#FFD700',       # Жёлтое предупреждение
}

class FastqReader:
    """
    FASTQ ридер для анализа файлов формата FASTQ
    """
    
    def __init__(self, filename):
        self.filename = filename
        self._sequence_count = None
        self._total_length = None
        self._quality_data = None
        self._length_data = None
        self._content_data = None
    
    def _read_fastq_chunks(self):
        """Генератор: читает FASTQ файл по одному риду за раз"""
        with open(self.filename, 'r', encoding='utf-8') as file:
            while True:
                lines = [file.readline().strip() for _ in range(4)]
                if not lines[0]:  # Конец файла
                    break
                yield lines
    
    def calculate_statistics(self):
        """Рассчитывает статистику используя генератор"""
        count = 0
        total_length = 0
        
        for chunk in self._read_fastq_chunks():
            count += 1
            sequence_line = chunk[1]  # Вторая строка - последовательность
            total_length += len(sequence_line)
        
        self._sequence_count = count
        self._total_length = total_length
        return count, total_length
    
    def get_sequence_count(self):
        """Возвращает количество последовательностей"""
        if self._sequence_count is None:
            self.calculate_statistics()
        return self._sequence_count
    
    def get_average_length(self):
        """Возвращает среднюю длину последовательностей"""
        if self._sequence_count is None:
            self.calculate_statistics()
        if self._sequence_count == 0:
            return 0
        return self._total_length / self._sequence_count
    
    def collect_quality_data(self):
        """Собирает данные для графика качества"""
        if self._quality_data is not None:
            return self._quality_data
            
        quality_sums = defaultdict(int)
        quality_counts = defaultdict(int)
        max_position = 0
        
        for chunk in self._read_fastq_chunks():
            quality_line = chunk[3]  # Четвертая строка - качество
            for i, char in enumerate(quality_line):
                score = ord(char) - 33  # Конвертация в числовое качество
                quality_sums[i] += score
                quality_counts[i] += 1
                max_position = max(max_position, i)
        
        positions = list(range(1, max_position + 2))
        avg_qualities = [quality_sums[i] / quality_counts[i] for i in range(max_position + 1)]
        
        self._quality_data = (positions, avg_qualities)
        return self._quality_data
    
    def collect_length_data(self):
        """Собирает данные для гистограммы длин"""
        if self._length_data is not None:
            return self._length_data
            
        lengths = []
        
        for chunk in self._read_fastq_chunks():
            sequence_line = chunk[1]  # Вторая строка - последовательность
            lengths.append(len(sequence_line))
        
        self._length_data = lengths
        return self._length_data
    
    def collect_content_data(self):
        """Собирает данные для графика содержания нуклеотидов"""
        if self._content_data is not None:
            return self._content_data
            
        base_counts = {'A': defaultdict(int), 'C': defaultdict(int), 
                      'G': defaultdict(int), 'T': defaultdict(int)}
        total_counts = defaultdict(int)
        max_position = 0
        
        for chunk in self._read_fastq_chunks():
            sequence_line = chunk[1]  # Вторая строка - последовательность
            for i, base in enumerate(sequence_line):
                base_upper = base.upper()
                if base_upper in base_counts:
                    base_counts[base_upper][i] += 1
                    total_counts[i] += 1
                    max_position = max(max_position, i)
        
        positions = list(range(1, max_position + 2))
        content_data = {}
        for base in base_counts:
            counts = base_counts[base]
            percentages = [counts[i] / total_counts[i] * 100 if total_counts[i] > 0 else 0 
                          for i in range(max_position + 1)]
            content_data[base] = percentages
        
        self._content_data = (positions, content_data)
        return self._content_data


class FastqAnalyzerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("🧬 FASTQ File Analyzer")
        self.root.geometry("1200x800")
        self.root.configure(bg=PISTACHIO_THEME['secondary'])
        
        self.current_file = None
        self.analyzer = None
        
        self.setup_styles()
        self.setup_gui()
        
    def setup_styles(self):
        """Настраиваем стили для виджетов"""
        self.style = ttk.Style()
        
        # Настраиваем тему
        self.style.theme_use('clam')
        
        # Стиль для кнопок
        self.style.configure('Pistachio.TButton',
                           background=PISTACHIO_THEME['primary'],
                           foreground=PISTACHIO_THEME['text_dark'],
                           borderwidth=0,
                           focuscolor=PISTACHIO_THEME['primary_light'],
                           padding=(15, 8))
        
        self.style.map('Pistachio.TButton',
                      background=[('active', PISTACHIO_THEME['primary_light']),
                                 ('pressed', PISTACHIO_THEME['primary_dark'])])
        
        # Стиль для фреймов
        self.style.configure('Pistachio.TLabelframe',
                           background=PISTACHIO_THEME['secondary'],
                           foreground=PISTACHIO_THEME['text_dark'],
                           borderwidth=2,
                           relief='groove')
        
        self.style.configure('Pistachio.TLabelframe.Label',
                           background=PISTACHIO_THEME['secondary'],
                           foreground=PISTACHIO_THEME['text_dark'],
                           font=('Arial', 10, 'bold'))
        
        # Стиль для прогресс-бара
        self.style.configure('Pistachio.Horizontal.TProgressbar',
                           background=PISTACHIO_THEME['primary'],
                           troughcolor=PISTACHIO_THEME['secondary'],
                           borderwidth=0,
                           lightcolor=PISTACHIO_THEME['primary_light'],
                           darkcolor=PISTACHIO_THEME['primary_dark'])
        
        # Стиль для вкладок
        self.style.configure('Pistachio.TNotebook',
                           background=PISTACHIO_THEME['secondary'],
                           borderwidth=0)
        
        self.style.configure('Pistachio.TNotebook.Tab',
                           background=PISTACHIO_THEME['primary_light'],
                           foreground=PISTACHIO_THEME['text_dark'],
                           padding=(20, 8),
                           font=('Arial', 9, 'bold'))
        
        self.style.map('Pistachio.TNotebook.Tab',
                      background=[('selected', PISTACHIO_THEME['primary']),
                                 ('active', PISTACHIO_THEME['primary_light'])])
    
    def setup_gui(self):
        # Header
        header_frame = tk.Frame(self.root, bg=PISTACHIO_THEME['primary'], height=80)
        header_frame.pack(fill=tk.X, padx=10, pady=(10, 5))
        header_frame.pack_propagate(False)
        
        header_label = tk.Label(header_frame, 
                               text="🧬 FASTQ File Analyzer 🧬", 
                               font=('Arial', 16, 'bold'),
                               bg=PISTACHIO_THEME['primary'],
                               fg=PISTACHIO_THEME['text_dark'])
        header_label.pack(expand=True)
        
        # Main frame
        main_frame = tk.Frame(self.root, bg=PISTACHIO_THEME['secondary'])
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # File selection section
        file_frame = ttk.LabelFrame(main_frame, text=" File Selection", style='Pistachio.TLabelframe')
        file_frame.pack(fill=tk.X, pady=(0, 10))
        
        file_content = tk.Frame(file_frame, bg=PISTACHIO_THEME['secondary'])
        file_content.pack(fill=tk.X, padx=10, pady=10)
        
        browse_btn = ttk.Button(file_content, text=" Browse FASTQ File", 
                               command=self.browse_file, style='Pistachio.TButton')
        browse_btn.pack(side=tk.LEFT, padx=(0, 10))
        
        self.file_label = tk.Label(file_content, text="No file selected", 
                                  font=('Arial', 9),
                                  bg=PISTACHIO_THEME['secondary'],
                                  fg=PISTACHIO_THEME['text_dark'],
                                  wraplength=800)
        self.file_label.pack(side=tk.LEFT, fill=tk.X, expand=True)
        
        # Analysis controls
        control_frame = ttk.LabelFrame(main_frame, text="🔧 Analysis Controls", style='Pistachio.TLabelframe')
        control_frame.pack(fill=tk.X, pady=(0, 10))
        
        control_content = tk.Frame(control_frame, bg=PISTACHIO_THEME['secondary'])
        control_content.pack(fill=tk.X, padx=10, pady=10)
        
        # Кнопки анализа
        buttons_frame = tk.Frame(control_content, bg=PISTACHIO_THEME['secondary'])
        buttons_frame.pack(fill=tk.X)
        
        ttk.Button(buttons_frame, text=" Show Statistics", 
                  command=self.show_statistics, style='Pistachio.TButton').pack(side=tk.LEFT, padx=(0, 8))
        
        ttk.Button(buttons_frame, text=" Quality Analysis", 
                  command=self.quality_analysis, style='Pistachio.TButton').pack(side=tk.LEFT, padx=(0, 8))
        
        ttk.Button(buttons_frame, text=" Nucleotide Content", 
                  command=self.nucleotide_analysis, style='Pistachio.TButton').pack(side=tk.LEFT, padx=(0, 8))
        
        ttk.Button(buttons_frame, text=" Full Analysis", 
                  command=self.full_analysis, style='Pistachio.TButton').pack(side=tk.LEFT)
        
        # Progress bar
        self.progress = ttk.Progressbar(control_content, mode='indeterminate', 
                                       style='Pistachio.Horizontal.TProgressbar')
        self.progress.pack(fill=tk.X, pady=(10, 0))
        
        # Results section
        results_frame = ttk.LabelFrame(main_frame, text=" Analysis Results", style='Pistachio.TLabelframe')
        results_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))
        
        # Create notebook for tabs
        self.notebook = ttk.Notebook(results_frame, style='Pistachio.TNotebook')
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Statistics tab
        self.stats_frame = tk.Frame(self.notebook, bg=PISTACHIO_THEME['secondary'])
        self.notebook.add(self.stats_frame, text=" Statistics")
        
        self.stats_text = tk.Text(self.stats_frame, wrap=tk.WORD, width=80, height=15,
                                 font=('Consolas', 10),
                                 bg=PISTACHIO_THEME['secondary'],
                                 fg=PISTACHIO_THEME['text_dark'],
                                 relief='flat',
                                 padx=10, pady=10)
        stats_scrollbar = ttk.Scrollbar(self.stats_frame, orient="vertical", command=self.stats_text.yview)
        self.stats_text.configure(yscrollcommand=stats_scrollbar.set)
        
        self.stats_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        stats_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Plots tab
        self.plots_frame = tk.Frame(self.notebook, bg=PISTACHIO_THEME['secondary'])
        self.notebook.add(self.plots_frame, text=" Plots")
        
        # Status bar
        status_frame = tk.Frame(self.root, bg=PISTACHIO_THEME['primary_dark'], height=25)
        status_frame.pack(fill=tk.X, padx=10, pady=(0, 10))
        status_frame.pack_propagate(False)
        
        self.status_var = tk.StringVar(value="🚀 Ready to analyze FASTQ files")
        status_bar = tk.Label(status_frame, textvariable=self.status_var, 
                             font=('Arial', 8),
                             bg=PISTACHIO_THEME['primary_dark'],
                             fg=PISTACHIO_THEME['text_light'])
        status_bar.pack(side=tk.LEFT, padx=10)
        
        # Footer with version
        version_label = tk.Label(status_frame, text="v1.0 • FASTQ Analyzer", 
                                font=('Arial', 7),
                                bg=PISTACHIO_THEME['primary_dark'],
                                fg=PISTACHIO_THEME['text_light'])
        version_label.pack(side=tk.RIGHT, padx=10)
    
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
    
    def load_file(self, filename):
        """Load and validate FASTQ file"""
        self.current_file = filename
        self.file_label.config(text=f" {os.path.basename(filename)}")
        self.status_var.set(f" Loaded: {os.path.basename(filename)}")
        
        # Clear previous results
        self.stats_text.delete(1.0, tk.END)
        self.clear_plots()
        self.analyzer = None
    
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
        stats_text = f"""🧬 FASTQ FILE STATISTICS
{'=' * 50}

 File: {os.path.basename(self.current_file)}

 Basic Statistics:
• Sequence count: {count:,}
• Average length: {avg_len:.2f} bp
• Total data volume: {total_bp:,.0f} bp

 File Information:
• File size: {os.path.getsize(self.current_file) / 1024 / 1024:.2f} MB
• Memory usage: Optimized (generators)

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
                
                # Собираем данные для графиков
                self.analyzer.collect_quality_data()
                self.analyzer.collect_length_data()
                
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
        
        if not self.analyzer:
            return
        
        # Create frame for plots
        plots_container = tk.Frame(self.plots_frame, bg=PISTACHIO_THEME['secondary'])
        plots_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        try:
            # График 1: Качество по позициям
            positions, avg_qualities = self.analyzer.collect_quality_data()
            
            fig1, ax1 = plt.subplots(figsize=(6, 4))
            ax1.plot(positions, avg_qualities, linewidth=2.5, color=PISTACHIO_THEME['primary_dark'], 
                    marker='o', markersize=2, alpha=0.8)
            ax1.fill_between(positions, avg_qualities, alpha=0.3, color=PISTACHIO_THEME['primary_light'])
            ax1.set_title('Per Base Sequence Quality', fontsize=11, fontweight='bold', 
                         color=PISTACHIO_THEME['text_dark'], pad=10)
            ax1.set_xlabel('Position (bp)', fontsize=9, color=PISTACHIO_THEME['text_dark'])
            ax1.set_ylabel('Quality Score', fontsize=9, color=PISTACHIO_THEME['text_dark'])
            ax1.grid(True, alpha=0.3)
            ax1.set_facecolor(PISTACHIO_THEME['secondary'])
            fig1.patch.set_facecolor(PISTACHIO_THEME['secondary'])
            fig1.tight_layout()
            
            canvas1 = FigureCanvasTkAgg(fig1, plots_container)
            canvas1.draw()
            canvas1.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)
            
            # График 2: Распределение длин
            lengths = self.analyzer.collect_length_data()
            
            fig2, ax2 = plt.subplots(figsize=(6, 4))
            ax2.hist(lengths, bins=15, edgecolor=PISTACHIO_THEME['primary_dark'], 
                    alpha=0.7, color=PISTACHIO_THEME['primary'], linewidth=1.2)
            ax2.set_title('Sequence Length Distribution', fontsize=11, fontweight='bold', 
                         color=PISTACHIO_THEME['text_dark'], pad=10)
            ax2.set_xlabel('Length (bp)', fontsize=9, color=PISTACHIO_THEME['text_dark'])
            ax2.set_ylabel('Frequency', fontsize=9, color=PISTACHIO_THEME['text_dark'])
            ax2.grid(True, alpha=0.3)
            ax2.set_facecolor(PISTACHIO_THEME['secondary'])
            fig2.patch.set_facecolor(PISTACHIO_THEME['secondary'])
            fig2.tight_layout()
            
            canvas2 = FigureCanvasTkAgg(fig2, plots_container)
            canvas2.draw()
            canvas2.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
            
            # Switch to plots tab
            self.notebook.select(1)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to create plots: {str(e)}")
    
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
                
                # Собираем данные для графика содержания
                self.analyzer.collect_content_data()
                
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
        
        if not self.analyzer:
            return
        
        try:
            # Собираем данные для графика содержания
            positions, content_data = self.analyzer.collect_content_data()
            
            # Создаём график
            fig, ax = plt.subplots(figsize=(10, 5))
            
            colors = [PISTACHIO_THEME['primary'], '#FF6B6B', '#4ECDC4', '#FFD166']
            bases = ['A', 'C', 'G', 'T']
            
            for base, color in zip(bases, colors):
                percentages = content_data[base]
                ax.plot(positions, percentages, label=base, linewidth=2, color=color, alpha=0.8)
            
            ax.set_title('Per Base Sequence Content', fontsize=12, fontweight='bold', 
                        color=PISTACHIO_THEME['text_dark'], pad=15)
            ax.set_xlabel('Position in read (bp)', fontsize=10, color=PISTACHIO_THEME['text_dark'])
            ax.set_ylabel('Percentage (%)', fontsize=10, color=PISTACHIO_THEME['text_dark'])
            ax.legend(frameon=True, facecolor=PISTACHIO_THEME['secondary'], fontsize=9)
            ax.grid(True, alpha=0.3)
            ax.set_facecolor(PISTACHIO_THEME['secondary'])
            fig.patch.set_facecolor(PISTACHIO_THEME['secondary'])
            fig.tight_layout()
            
            canvas = FigureCanvasTkAgg(fig, self.plots_frame)
            canvas.draw()
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
            
            self.notebook.select(1)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to create nucleotide plot: {str(e)}")
    
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
                
                # Собираем данные для всех графиков
                self.analyzer.collect_quality_data()
                self.analyzer.collect_length_data()
                self.analyzer.collect_content_data()
                
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
        
        messagebox.showinfo("Analysis Complete", 
                          " Full analysis completed successfully!\n\n"
                          "Check the Statistics and Plots tabs for detailed results.")
    
    def start_processing(self):
        """Start progress indicator"""
        self.progress.start()
        self.status_var.set(" Processing... Please wait")
    
    def stop_processing(self):
        """Stop progress indicator"""
        self.progress.stop()
        self.status_var.set(" Analysis complete")

def main():
    # Create and run the GUI application
    root = tk.Tk()
    app = FastqAnalyzerGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
    root.mainloop()

if __name__ == "__main__":
    main()
