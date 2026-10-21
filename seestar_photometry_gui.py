"""Windows GUI for the unchanged Seestar photometry processing engine."""
from __future__ import annotations
import contextlib
import io
import os
import queue
import sys
import threading
import traceback
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import seestar_photometry as engine


class QueueWriter(io.TextIOBase):
    def __init__(self, messages):
        self.messages = messages
    def write(self, text):
        if text:
            self.messages.put(('log', text))
        return len(text)
    def flush(self):
        pass


class Application(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('Seestar Photometry')
        self.geometry('780x580')
        self.minsize(680, 480)
        self.messages = queue.Queue()
        self.folder = tk.StringVar()
        self.star = tk.StringVar()
        self.radius = tk.StringVar(value='8.0')
        self.running = False
        self.build_ui()
        self.after(100, self.drain_messages)

    def build_ui(self):
        outer = ttk.Frame(self, padding=18)
        outer.pack(fill='both', expand=True)
        ttk.Label(outer, text='Seestar Photometry', font=('Segoe UI', 18, 'bold')).pack(anchor='w', pady=(0, 12))
        ttk.Label(outer, text='CSV-adatok könyvtára').pack(anchor='w')
        folder_row = ttk.Frame(outer)
        folder_row.pack(fill='x', pady=(5, 12))
        ttk.Entry(folder_row, textvariable=self.folder).pack(side='left', fill='x', expand=True)
        ttk.Button(folder_row, text='Tallózás...', command=self.pick_folder).pack(side='left', padx=(8, 0))
        ttk.Button(folder_row, text='Frissítés', command=self.refresh_stars).pack(side='left', padx=(8, 0))
        options = ttk.Frame(outer)
        options.pack(fill='x', pady=(0, 12))
        ttk.Label(options, text='Célcsillag:').pack(side='left')
        self.star_select = ttk.Combobox(options, textvariable=self.star, state='readonly', width=32)
        self.star_select.pack(side='left', padx=(8, 20))
        ttk.Label(options, text='Keresési sugár (″):').pack(side='left')
        ttk.Entry(options, textvariable=self.radius, width=8).pack(side='left', padx=(8, 0))
        ttk.Label(outer, text='Csak a kijelölt könyvtár CSV-fájljait dolgozza fel. Az eredményeket az Eredmenyek almappába menti.', wraplength=720).pack(anchor='w', pady=(0, 12))
        action_row = ttk.Frame(outer)
        action_row.pack(fill='x', pady=(0, 14))
        self.start_btn = ttk.Button(action_row, text='Feldolgozás indítása', command=self.start)
        self.start_btn.pack(side='left')
        ttk.Button(action_row, text='Célkönyvtár megnyitása', command=self.open_results_folder).pack(side='left', padx=(10, 0))
        ttk.Label(outer, text='Feldolgozási napló').pack(anchor='w', pady=(0, 5))
        log_frame = ttk.Frame(outer)
        log_frame.pack(fill='both', expand=True)
        self.log = tk.Text(log_frame, wrap='word', state='disabled', height=12, font=('Consolas', 9))
        self.log.pack(side='left', fill='both', expand=True)
        scroll = ttk.Scrollbar(log_frame, command=self.log.yview)
        scroll.pack(side='right', fill='y')
        self.log.configure(yscrollcommand=scroll.set)

    def pick_folder(self):
        folder = filedialog.askdirectory(title='Válaszd ki a Seestar CSV könyvtárát')
        if folder:
            self.folder.set(folder)
            self.refresh_stars()

    def refresh_stars(self):
        folder = Path(self.folder.get().strip())
        if not folder.is_dir():
            messagebox.showwarning('Könyvtár', 'Válassz létező könyvtárat.')
            return
        stars = []
        for path in folder.glob('*.txt'):
            if path.name.lower() == 'readme.txt':
                continue
            try:
                cfg = engine.read_targets(path)
                if len(cfg.sections()) == 1 and path.stem.strip().casefold() == cfg.sections()[0].strip().casefold():
                    stars.append(cfg.sections()[0])
            except Exception:
                continue
        stars.sort(key=str.casefold)
        self.star_select['values'] = stars
        self.star.set(stars[0] if stars else '')
        if not stars:
            messagebox.showwarning('Csillagadatlap', 'Nem található érvényes csillagnév.txt a könyvtárban.')

    def open_results_folder(self):
        """Open the output folder beside the selected CSV files in Windows Explorer."""
        raw_folder = self.folder.get().strip()
        if not raw_folder:
            messagebox.showwarning('Könyvtár', 'Előbb válaszd ki a CSV-adatok könyvtárát.')
            return
        folder = Path(raw_folder)
        if not folder.is_dir():
            messagebox.showwarning('Könyvtár', 'A kiválasztott könyvtár nem található.')
            return
        output = folder / 'Eredmenyek'
        if not output.is_dir():
            messagebox.showinfo('Nincs még eredmény',
                                'Az Eredmenyek almappa még nem létezik. '
                                'A feldolgozás után automatikusan létrejön.')
            return
        try:
            os.startfile(str(output))
        except OSError as exc:
            messagebox.showerror('Megnyitási hiba', str(exc))

    def append(self, message):
        self.log.configure(state='normal')
        self.log.insert('end', message)
        self.log.see('end')
        self.log.configure(state='disabled')

    def start(self):
        folder = Path(self.folder.get().strip())
        if not folder.is_dir():
            messagebox.showerror('Hiba', 'A kiválasztott könyvtár nem létezik.')
            return
        star = self.star.get().strip()
        if not star:
            messagebox.showerror('Hiba', 'Válassz célcsillagot, szükség esetén a Frissítés gombbal.')
            return
        try:
            radius = float(self.radius.get().strip().replace(',', '.'))
            if not 0 < radius <= 3600:
                raise ValueError
        except ValueError:
            messagebox.showerror('Hiba', 'A keresési sugár pozitív szám legyen.')
            return
        self.running = True
        self.start_btn.configure(state='disabled')
        self.append(f'\n=== {star} | {folder} ===\n')
        threading.Thread(target=self.process, args=(folder, star, radius), daemon=True).start()

    def process(self, folder, star, radius):
        # The existing main() is preserved; only the command-line arguments are supplied by the GUI.
        previous = sys.argv
        writer = QueueWriter(self.messages)
        try:
            sys.argv = ['seestar_photometry.py', '--folder', str(folder), '--star', star, '--max-sep', str(radius)]
            with contextlib.redirect_stdout(writer), contextlib.redirect_stderr(writer):
                engine.main()
            self.messages.put(('done', True))
        except BaseException:
            self.messages.put(('log', '\nHIBA:\n' + traceback.format_exc() + '\n'))
            self.messages.put(('done', False))
        finally:
            sys.argv = previous

    def drain_messages(self):
        try:
            while True:
                action, value = self.messages.get_nowait()
                if action == 'log':
                    self.append(value)
                elif action == 'done':
                    self.running = False
                    self.start_btn.configure(state='normal')
                    if value:
                        self.append('\nFeldolgozás befejezve.\n')
                        messagebox.showinfo('Kész', 'A feldolgozás befejeződött. Az eredmények az Eredmenyek mappában vannak.')
                    else:
                        messagebox.showerror('Hiba', 'A feldolgozás során hiba történt. Részletek a naplóban.')
        except queue.Empty:
            pass
        self.after(100, self.drain_messages)


if __name__ == '__main__':
    Application().mainloop()
