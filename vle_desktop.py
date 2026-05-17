"""
VLE Desktop — Equilíbrio Líquido-Vapor Bicomponente
Interface CustomTkinter + Matplotlib
Banco de dados externo: componentes.json e nrtl.json
"""

import math, json, sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import customtkinter as ctk
from tkinter import filedialog, messagebox
import tkinter as tk
import tkinter.ttk as ttk

# ─────────────────────────────────────────────
# Localiza JSONs na mesma pasta do exe/script
# ─────────────────────────────────────────────
def base_dir():
    if getattr(sys, "frozen", False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def carregar_componentes():
    path = os.path.join(base_dir(), "componentes.json")
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except FileNotFoundError:
        messagebox.showerror("Erro",
            f"Arquivo não encontrado:\n{path}\n\n"
            "Coloque 'componentes.json' na mesma pasta do programa.")
        sys.exit(1)

def carregar_nrtl():
    path = os.path.join(base_dir(), "nrtl.json")
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except FileNotFoundError:
        messagebox.showerror("Erro",
            f"Arquivo não encontrado:\n{path}\n\n"
            "Coloque 'nrtl.json' na mesma pasta do programa.")
        sys.exit(1)

COMPONENTS = carregar_componentes()
NRTL_RAW   = carregar_nrtl()
COMP_NAMES = sorted(COMPONENTS.keys())

def get_nrtl_params(name1, name2):
    key  = f"{name1}|{name2}"
    keyT = f"{name2}|{name1}"
    if key in NRTL_RAW:
        p = NRTL_RAW[key]
        return p["A12"], p["A21"], p["alpha"]
    elif keyT in NRTL_RAW:
        p = NRTL_RAW[keyT]
        return p["A21"], p["A12"], p["alpha"]
    return None, None, None

# ─────────────────────────────────────────────
# Tema
# ─────────────────────────────────────────────
ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

# ─────────────────────────────────────────────
# Termodinâmica
# ─────────────────────────────────────────────
def pvap(A, B, C, T_K):
    return 10.0 ** (A - B / (T_K + C))

def nrtl_gamma(x1, T_K, A12, A21, alpha):
    R  = 8.314
    x2 = 1.0 - x1
    tau12 = A12 / (R * T_K)
    tau21 = A21 / (R * T_K)
    G12   = math.exp(-alpha * tau12)
    G21   = math.exp(-alpha * tau21)
    d1 = x1 + x2 * G21
    d2 = x2 + x1 * G12
    eps = 1e-12
    if abs(d1) < eps: d1 = eps
    if abs(d2) < eps: d2 = eps
    lng1 = x2**2 * (tau21*(G21/d1)**2 + tau12*G12/d2**2)
    lng2 = x1**2 * (tau12*(G12/d2)**2 + tau21*G21/d1**2)
    return math.exp(lng1), math.exp(lng2)

def bubble_T(x1, A1, B1, C1, A2, B2, C2, P_bar, T_init=None,
             modelo="Raoult", A12=0.0, A21=0.0, alpha_nrtl=0.30):
    if T_init is None:
        Tb1 = B1 / (A1 - math.log10(P_bar)) - C1
        Tb2 = B2 / (A2 - math.log10(P_bar)) - C2
        T   = x1*Tb1 + (1-x1)*Tb2
    else:
        T = T_init
    for _ in range(200):
        pb1 = pvap(A1, B1, C1, T)
        pb2 = pvap(A2, B2, C2, T)
        g1, g2 = (nrtl_gamma(x1, T, A12, A21, alpha_nrtl)
                  if modelo == "NRTL" else (1.0, 1.0))
        Pcalc = x1*g1*pb1 + (1-x1)*g2*pb2
        err   = Pcalc - P_bar
        if abs(err) < 1e-7: break
        dT   = 0.01
        pb1d = pvap(A1, B1, C1, T+dT)
        pb2d = pvap(A2, B2, C2, T+dT)
        g1d, g2d = (nrtl_gamma(x1, T+dT, A12, A21, alpha_nrtl)
                    if modelo == "NRTL" else (1.0, 1.0))
        dPdT = (x1*g1d*pb1d + (1-x1)*g2d*pb2d - Pcalc) / dT
        if abs(dPdT) < 1e-12: break
        T -= err / dPdT
    return T

def calc_vle(A1, B1, C1, A2, B2, C2, P_bar, n_points,
             modelo="Raoult", A12=0.0, A21=0.0, alpha_nrtl=0.30):
    rows, T_prev = [], None
    for x1 in np.linspace(0, 1, n_points):
        T = bubble_T(x1, A1, B1, C1, A2, B2, C2, P_bar,
                     T_init=T_prev, modelo=modelo,
                     A12=A12, A21=A21, alpha_nrtl=alpha_nrtl)
        T_prev = T
        P1s = pvap(A1, B1, C1, T)
        P2s = pvap(A2, B2, C2, T)
        g1, g2 = (nrtl_gamma(x1, T, A12, A21, alpha_nrtl)
                  if modelo == "NRTL" else (1.0, 1.0))
        y1    = max(0.0, min(1.0, x1*g1*P1s/P_bar))
        alpha = (g1*P1s)/(g2*P2s) if P2s > 0 and g2 > 0 else np.nan
        rows.append({
            "x₁": round(float(x1),6),   "y₁": round(float(y1),6),
            "T (K)": round(float(T),4), "T (°C)": round(float(T-273.15),4),
            "P₁ˢᵃᵗ (bar)": round(float(P1s),6),
            "P₂ˢᵃᵗ (bar)": round(float(P2s),6),
            "γ₁": round(float(g1),6),   "γ₂": round(float(g2),6),
            "α₁₂ efetivo": round(float(alpha),4),
        })
    return pd.DataFrame(rows)

def detectar_azeotropo(df):
    diff = df["y₁"].values - df["x₁"].values
    for i in range(len(diff)-1):
        if diff[i]*diff[i+1] < 0:
            x0, x1_ = df["x₁"].iloc[i], df["x₁"].iloc[i+1]
            d0, d1  = diff[i], diff[i+1]
            frac    = -d0/(d1-d0)
            x_az    = x0 + frac*(x1_-x0)
            T_az    = (df["T (°C)"].iloc[i]
                       + frac*(df["T (°C)"].iloc[i+1]-df["T (°C)"].iloc[i]))
            return x_az, x_az, T_az
    return None

# ─────────────────────────────────────────────
# Cores dos gráficos
# ─────────────────────────────────────────────
AXBG="#1A1D23"; GRID="#2A2D35"; TXT="#E0DDD6"
BLUE="#4A9EDF"; ORANGE="#E87040"; GREEN="#5BB87A"; GRAY="#888880"

def style_ax(ax, title, xlabel, ylabel):
    ax.set_facecolor(AXBG)
    ax.tick_params(colors=TXT, labelsize=9)
    ax.xaxis.label.set_color(TXT); ax.yaxis.label.set_color(TXT)
    ax.title.set_color(TXT)
    ax.set_title(title, fontsize=10, fontweight="bold", pad=6)
    ax.set_xlabel(xlabel, fontsize=9); ax.set_ylabel(ylabel, fontsize=9)
    ax.grid(color=GRID, linewidth=0.6, linestyle="--")
    for sp in ax.spines.values(): sp.set_edgecolor(GRID)

# ─────────────────────────────────────────────
# Widget de autocompletar
# ─────────────────────────────────────────────
class AutocompleteEntry(ctk.CTkFrame):
    def __init__(self, master, names, on_select, width=270, **kwargs):
        super().__init__(master, fg_color="transparent", **kwargs)
        self._names     = names
        self._on_select = on_select
        self._popup     = None

        self.entry = ctk.CTkEntry(self, width=width,
                                   placeholder_text="Digite para buscar...")
        self.entry.pack()
        self.entry.bind("<KeyRelease>", self._on_key)
        self.entry.bind("<FocusOut>",   lambda e: self.after(150, self._fechar_popup))
        self.entry.bind("<Escape>",     lambda e: self._fechar_popup())

    def get(self):
        return self.entry.get()

    def set(self, v):
        self.entry.delete(0, "end")
        self.entry.insert(0, v)

    def _on_key(self, event):
        typed = self.entry.get().strip().lower()
        if not typed:
            self._fechar_popup(); return
        sugestoes = [n for n in self._names if typed in n.lower()]
        if sugestoes: self._mostrar_popup(sugestoes)
        else:         self._fechar_popup()

    def _mostrar_popup(self, sugestoes):
        self._fechar_popup()
        popup = tk.Toplevel(self)
        popup.wm_overrideredirect(True)
        popup.configure(bg="#1A1D23")
        self.entry.update_idletasks()
        x = self.entry.winfo_rootx()
        y = self.entry.winfo_rooty() + self.entry.winfo_height()
        w = self.entry.winfo_width()
        h = min(len(sugestoes)*28, 220)
        popup.geometry(f"{w}x{h}+{x}+{y}")

        frame = tk.Frame(popup, bg="#1A1D23")
        frame.pack(fill="both", expand=True)
        sb = tk.Scrollbar(frame); sb.pack(side="right", fill="y")
        lb = tk.Listbox(frame, bg="#1A1D23", fg="#E0DDD6",
                        selectbackground="#4A9EDF", selectforeground="white",
                        font=("Consolas", 11), bd=0, highlightthickness=0,
                        yscrollcommand=sb.set, activestyle="none")
        lb.pack(side="left", fill="both", expand=True)
        sb.config(command=lb.yview)
        for s in sugestoes: lb.insert("end", s)

        def selecionar(e=None):
            sel = lb.curselection()
            if sel:
                nome = lb.get(sel[0])
                self.set(nome)
                self._fechar_popup()
                self._on_select(nome)

        lb.bind("<ButtonRelease-1>", selecionar)
        lb.bind("<Return>",          selecionar)
        self._popup = popup

    def _fechar_popup(self):
        if self._popup:
            try: self._popup.destroy()
            except: pass
            self._popup = None

# ─────────────────────────────────────────────
# App principal
# ─────────────────────────────────────────────
class VLEApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("⚗️ Equilíbrio VLE Bicomponente")
        self.geometry("1400x860")
        self.minsize(1100, 700)
        self.configure(fg_color="#0E1117")

        self.df_result = None
        self.df_raoult = None
        self._c1_label = "Benzeno"
        self._c2_label = "Tolueno"
        self._c1_data  = COMPONENTS.get("Benzeno", {})
        self._c2_data  = COMPONENTS.get("Tolueno",  {})

        self._build_ui()
        self.after(200, self._calcular)

    # ── Helpers ───────────────────────────────
    def _section(self, parent, text):
        ctk.CTkLabel(parent, text=text, font=ctk.CTkFont(size=12, weight="bold"),
                     text_color="#4A9EDF").pack(anchor="w", padx=12, pady=(14,2))
        ctk.CTkFrame(parent, height=1, fg_color="#2A2D35").pack(fill="x", padx=8, pady=(0,6))

    def _lbl(self, parent, text):
        ctk.CTkLabel(parent, text=text, font=ctk.CTkFont(size=11),
                     text_color="#AAAAAA").pack(anchor="w", padx=14, pady=(4,0))

    # ── Build UI ──────────────────────────────
    def _build_ui(self):
        self.sidebar = ctk.CTkScrollableFrame(self, width=310, fg_color="#13161D",
                                               corner_radius=0)
        self.sidebar.pack(side="left", fill="y")
        self.main_area = ctk.CTkFrame(self, fg_color="#0E1117", corner_radius=0)
        self.main_area.pack(side="left", fill="both", expand=True)
        self._build_sidebar()
        self._build_main_area()

    def _build_sidebar(self):
        sb = self.sidebar

        ctk.CTkLabel(sb, text="⚗️  VLE  Bicomponente",
                     font=ctk.CTkFont(size=15, weight="bold"),
                     text_color="#E0DDD6").pack(pady=(16,4), padx=12, anchor="w")
        ctk.CTkLabel(sb, text="Raoult · NRTL · Antoine (DDB/NIST)",
                     font=ctk.CTkFont(size=10), text_color="#555555").pack(anchor="w", padx=14)

        # Modelo
        self._section(sb, "Modelo termodinâmico")
        self.modelo_var = ctk.StringVar(value="Raoult (ideal)")
        for opt in ["Raoult (ideal)", "NRTL (não ideal)"]:
            ctk.CTkRadioButton(sb, text=opt, variable=self.modelo_var, value=opt,
                               command=self._on_modelo_change,
                               font=ctk.CTkFont(size=11)).pack(anchor="w", padx=16, pady=2)

        # Comp 1
        self._section(sb, "Light Key (mais volátil)")
        self.ac1 = AutocompleteEntry(sb, COMP_NAMES, self._on_comp1_select)
        self.ac1.pack(padx=12, pady=4)
        self.ac1.set("Benzeno")
        self.c1_info = ctk.CTkLabel(sb, text="", font=ctk.CTkFont(size=10),
                                     text_color="#555555", wraplength=270)
        self.c1_info.pack(anchor="w", padx=14)
        self.c1_manual = ctk.CTkFrame(sb, fg_color="transparent")
        for lbl, attr in [("A₁","c1_A"),("B₁","c1_B"),("C₁","c1_C")]:
            self._lbl(self.c1_manual, lbl)
            e = ctk.CTkEntry(self.c1_manual, width=270); e.pack(padx=12, pady=2)
            setattr(self, attr, e)

        # Comp 2
        self._section(sb, "Heavy Key (menos volátil)")
        self.ac2 = AutocompleteEntry(sb, COMP_NAMES, self._on_comp2_select)
        self.ac2.pack(padx=12, pady=4)
        self.ac2.set("Tolueno")
        self.c2_info = ctk.CTkLabel(sb, text="", font=ctk.CTkFont(size=10),
                                     text_color="#555555", wraplength=270)
        self.c2_info.pack(anchor="w", padx=14)
        self.c2_manual = ctk.CTkFrame(sb, fg_color="transparent")
        for lbl, attr in [("A₂","c2_A"),("B₂","c2_B"),("C₂","c2_C")]:
            self._lbl(self.c2_manual, lbl)
            e = ctk.CTkEntry(self.c2_manual, width=270); e.pack(padx=12, pady=2)
            setattr(self, attr, e)

        # NRTL
        self.nrtl_frame = ctk.CTkFrame(sb, fg_color="transparent")
        self._section(self.nrtl_frame, "Parâmetros NRTL")
        self.nrtl_status = ctk.CTkLabel(self.nrtl_frame, text="",
                                         font=ctk.CTkFont(size=10),
                                         text_color="#5BB87A", wraplength=270)
        self.nrtl_status.pack(anchor="w", padx=14, pady=2)
        for lbl, attr, default in [("A₁₂ (J/mol)","nrtl_A12","1000.0"),
                                    ("A₂₁ (J/mol)","nrtl_A21","1000.0"),
                                    ("α (0.10–0.60)","nrtl_alpha","0.30")]:
            self._lbl(self.nrtl_frame, lbl)
            e = ctk.CTkEntry(self.nrtl_frame, width=270)
            e.insert(0, default); e.pack(padx=12, pady=2)
            setattr(self, attr, e)

        # Condições
        self._section(sb, "Condições de operação")
        self._lbl(sb, "Pressão (atm)")
        self.P_entry = ctk.CTkEntry(sb, width=270)
        self.P_entry.insert(0, "1.0"); self.P_entry.pack(padx=12, pady=2)

        self._lbl(sb, "Número de pontos")
        self.n_var = ctk.IntVar(value=21)
        self.n_slider = ctk.CTkSlider(sb, from_=11, to=101, number_of_steps=18,
                                       variable=self.n_var, width=250,
                                       command=lambda v: self.n_lbl.configure(
                                           text=f"{int(float(v))} pontos"))
        self.n_slider.pack(padx=16, pady=4)
        self.n_lbl = ctk.CTkLabel(sb, text="21 pontos", font=ctk.CTkFont(size=10),
                                   text_color="#AAAAAA")
        self.n_lbl.pack(anchor="w", padx=14)

        # Botões
        ctk.CTkButton(sb, text="🔄  Calcular", height=40,
                       font=ctk.CTkFont(size=13, weight="bold"),
                       fg_color="#4A9EDF", hover_color="#3A8ECF",
                       command=self._calcular).pack(padx=12, pady=(16,6), fill="x")

        for txt, cmd in [
            ("⬇  CSV completo",              self._export_csv_full),
            ("⬇  x,y para McCabe-Thiele",    self._export_csv_xy),
            ("🖼  Salvar gráficos (PNG)",     self._save_fig),
            ("🔄  Recarregar banco de dados", self._recarregar_banco),
        ]:
            ctk.CTkButton(sb, text=txt, height=34, font=ctk.CTkFont(size=11),
                           fg_color="#2A2D35", hover_color="#3A3D45",
                           command=cmd).pack(padx=12, pady=3, fill="x")

        ctk.CTkLabel(sb, text=f"📁 {base_dir()}",
                     font=ctk.CTkFont(size=9), text_color="#444444",
                     wraplength=270).pack(anchor="w", padx=14, pady=(8,16))

        self._on_comp1_select("Benzeno")
        self._on_comp2_select("Tolueno")
        self._on_modelo_change()

    def _build_main_area(self):
        ma = self.main_area

        # Métricas
        self.mf = ctk.CTkFrame(ma, height=60, fg_color="#13161D", corner_radius=8)
        self.mf.pack(fill="x", padx=12, pady=(10,0))
        self.mf.pack_propagate(False)
        self.mlbls = {}
        for key in ["P operação","T_eb comp1","T_eb comp2","α₁₂ médio","Azeótropo"]:
            f = ctk.CTkFrame(self.mf, fg_color="#1A1D23", corner_radius=6)
            f.pack(side="left", expand=True, fill="both", padx=6, pady=8)
            ctk.CTkLabel(f, text=key, font=ctk.CTkFont(size=9),
                         text_color="#555555").pack(pady=(4,0))
            l = ctk.CTkLabel(f, text="—", font=ctk.CTkFont(size=12, weight="bold"),
                              text_color="#E0DDD6")
            l.pack()
            self.mlbls[key] = l

        # Abas
        self.tabview = ctk.CTkTabview(ma, fg_color="#0E1117",
                                       segmented_button_fg_color="#13161D",
                                       segmented_button_selected_color="#4A9EDF",
                                       segmented_button_selected_hover_color="#3A8ECF",
                                       segmented_button_unselected_color="#13161D")
        self.tabview.pack(fill="both", expand=True, padx=12, pady=8)
        tab_g = self.tabview.add("📈  Gráficos")
        tab_t = self.tabview.add("📋  Tabela")

        self.fig = plt.figure(figsize=(13,8), facecolor="#0E1117")
        self.canvas = FigureCanvasTkAgg(self.fig, master=tab_g)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self._build_table(tab_t)

    def _build_table(self, parent):
        frame = ctk.CTkFrame(parent, fg_color="#0E1117")
        frame.pack(fill="both", expand=True, padx=8, pady=8)
        cols = ["x₁","y₁","T (K)","T (°C)","P₁ˢᵃᵗ (bar)","P₂ˢᵃᵗ (bar)","γ₁","γ₂","α₁₂ efetivo"]
        style = ttk.Style()
        style.theme_use("default")
        style.configure("Treeview", background="#1A1D23", foreground="#E0DDD6",
                         fieldbackground="#1A1D23", rowheight=24, font=("Consolas",10))
        style.configure("Treeview.Heading", background="#2A2D35", foreground="#4A9EDF",
                         font=("Consolas",10,"bold"), relief="flat")
        style.map("Treeview", background=[("selected","#4A9EDF")])
        self.tree = ttk.Treeview(frame, columns=cols, show="headings", height=30)
        for c in cols:
            self.tree.heading(c, text=c)
            self.tree.column(c, width=110, anchor="center")
        vsb = ttk.Scrollbar(frame, orient="vertical",   command=self.tree.yview)
        hsb = ttk.Scrollbar(frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)

    # ── Callbacks ─────────────────────────────
    def _on_comp1_select(self, nome):
        d = COMPONENTS.get(nome)
        if d:
            self._c1_label = nome
            self._c1_data  = d
            self.c1_info.configure(
                text=f"A={d['A']}  B={d['B']}  C={d['C']}  |  Tb={d['Tb_C']} °C")
            self.c1_manual.pack_forget()
        else:
            self._c1_label = nome or "Comp A"
            self._c1_data  = {}
            self.c1_info.configure(text="Componente não encontrado — preencha abaixo")
            self.c1_manual.pack(fill="x")
        self._atualizar_nrtl_status()

    def _on_comp2_select(self, nome):
        d = COMPONENTS.get(nome)
        if d:
            self._c2_label = nome
            self._c2_data  = d
            self.c2_info.configure(
                text=f"A={d['A']}  B={d['B']}  C={d['C']}  |  Tb={d['Tb_C']} °C")
            self.c2_manual.pack_forget()
        else:
            self._c2_label = nome or "Comp B"
            self._c2_data  = {}
            self.c2_info.configure(text="Componente não encontrado — preencha abaixo")
            self.c2_manual.pack(fill="x")
        self._atualizar_nrtl_status()

    def _on_modelo_change(self):
        if self.modelo_var.get() == "NRTL (não ideal)":
            self.nrtl_frame.pack(fill="x")
            self._atualizar_nrtl_status()
        else:
            self.nrtl_frame.pack_forget()

    def _atualizar_nrtl_status(self):
        A12, A21, alpha = get_nrtl_params(self._c1_label, self._c2_label)
        if A12 is not None:
            self.nrtl_status.configure(
                text=f"✅ Parâmetros encontrados: {self._c1_label} / {self._c2_label}",
                text_color="#5BB87A")
            self.nrtl_A12.delete(0,"end");   self.nrtl_A12.insert(0, str(A12))
            self.nrtl_A21.delete(0,"end");   self.nrtl_A21.insert(0, str(A21))
            self.nrtl_alpha.delete(0,"end"); self.nrtl_alpha.insert(0, str(alpha))
        else:
            self.nrtl_status.configure(
                text="⚠️ Par não encontrado no nrtl.json — insira manualmente",
                text_color="#E87040")

    # ── Parâmetros ────────────────────────────
    def _get_params(self):
        d1 = self._c1_data
        d2 = self._c2_data
        A1,B1,C1 = ((d1["A"],d1["B"],d1["C"]) if d1
                    else (float(self.c1_A.get()),float(self.c1_B.get()),float(self.c1_C.get())))
        A2,B2,C2 = ((d2["A"],d2["B"],d2["C"]) if d2
                    else (float(self.c2_A.get()),float(self.c2_B.get()),float(self.c2_C.get())))
        P_atm   = float(self.P_entry.get())
        n_pts   = int(self.n_var.get())
        nrtl    = self.modelo_var.get() == "NRTL (não ideal)"
        A12_v   = float(self.nrtl_A12.get())   if nrtl else 0.0
        A21_v   = float(self.nrtl_A21.get())   if nrtl else 0.0
        alpha_v = float(self.nrtl_alpha.get()) if nrtl else 0.30
        return (self._c1_label, A1,B1,C1,
                self._c2_label, A2,B2,C2,
                P_atm, n_pts, nrtl, A12_v, A21_v, alpha_v)

    # ── Cálculo ───────────────────────────────
    def _calcular(self):
        try:
            (c1,A1,B1,C1, c2,A2,B2,C2,
             P_atm,n_pts,nrtl,A12_v,A21_v,alpha_v) = self._get_params()

            P_bar = P_atm
            Tb1 = bubble_T(1.0, A1,B1,C1, A2,B2,C2, P_bar)
            Tb2 = bubble_T(0.0, A1,B1,C1, A2,B2,C2, P_bar)
            swapped = False
            if Tb1 > Tb2:
                A1,B1,C1, A2,B2,C2 = A2,B2,C2, A1,B1,C1
                c1,c2 = c2,c1
                A12_v,A21_v = A21_v,A12_v
                swapped = True

            mod = "NRTL" if nrtl else "Raoult"
            df  = calc_vle(A1,B1,C1, A2,B2,C2, P_bar, n_pts,
                           modelo=mod, A12=A12_v, A21=A21_v, alpha_nrtl=alpha_v)
            df_r = (calc_vle(A1,B1,C1, A2,B2,C2, P_bar, n_pts, modelo="Raoult")
                    if nrtl else None)
            az = detectar_azeotropo(df)
            self.df_result = df
            self.df_raoult = df_r

            am   = df["α₁₂ efetivo"].mean()
            Teb1 = df["T (°C)"].iloc[-1]
            Teb2 = df["T (°C)"].iloc[0]
            self.mlbls["P operação"].configure(text=f"{P_atm:.2f} atm")
            self.mlbls["T_eb comp1"].configure(text=f"{Teb1:.1f} °C  ({c1})")
            self.mlbls["T_eb comp2"].configure(text=f"{Teb2:.1f} °C  ({c2})")
            self.mlbls["α₁₂ médio"].configure(text=f"{am:.3f}")
            if az:
                self.mlbls["Azeótropo"].configure(
                    text=f"x≈{az[0]:.3f}  T≈{az[2]:.1f}°C", text_color="#E87040")
            else:
                self.mlbls["Azeótropo"].configure(text="Não detectado", text_color="#5BB87A")

            self._plotar(df, df_r, c1, c2, P_atm, nrtl, az)
            self._preencher_tabela(df)

            if swapped:
                messagebox.showinfo("Reordenação automática",
                    f"Componentes reordenados:\nLight Key = {c1}\nHeavy Key = {c2}")
        except Exception as e:
            messagebox.showerror("Erro no cálculo", str(e))

    # ── Gráficos ──────────────────────────────
    def _plotar(self, df, df_r, c1, c2, P_atm, nrtl, az):
        self.fig.clear()
        self.fig.patch.set_facecolor("#0E1117")
        gs = gridspec.GridSpec(2,2, figure=self.fig, hspace=0.40, wspace=0.32,
                                left=0.07, right=0.97, top=0.93, bottom=0.07)
        x=df["x₁"].values; y=df["y₁"].values; T=df["T (°C)"].values
        P1=df["P₁ˢᵃᵗ (bar)"].values; P2=df["P₂ˢᵃᵗ (bar)"].values
        al=df["α₁₂ efetivo"].values; g1=df["γ₁"].values; g2=df["γ₂"].values
        lm  = "NRTL" if nrtl else "Raoult"
        COL = ORANGE if nrtl else BLUE

        # y–x
        ax1 = self.fig.add_subplot(gs[0,0])
        style_ax(ax1, f"Diagrama y–x  ({c1} / {c2})", f"x₁  [{c1}]", f"y₁  [{c1}]")
        if nrtl and df_r is not None:
            ax1.plot(df_r["x₁"].values, df_r["y₁"].values,
                     color=BLUE, lw=1.5, ls="--", label="Raoult (ref.)", alpha=0.6)
        ax1.plot(x, y, color=COL, lw=2, label=lm)
        ax1.plot([0,1],[0,1], color=GRAY, lw=1, ls="--", label="y=x")
        ax1.scatter(x[1:-1], y[1:-1], color=COL, s=18, zorder=5)
        if az:
            ax1.scatter([az[0]],[az[1]], color="red", s=80, zorder=10,
                        marker="*", label=f"Azeótropo x≈{az[0]:.3f}")
        ax1.set_xlim(0,1); ax1.set_ylim(0,1)
        ax1.legend(fontsize=8, facecolor=AXBG, labelcolor=TXT, edgecolor=GRID)

        # T–x–y
        ax2 = self.fig.add_subplot(gs[0,1])
        style_ax(ax2, f"Diagrama T–x–y  ({P_atm:.2f} atm)", f"Fração molar [{c1}]", "T (°C)")
        ax2.plot(x, T, color=ORANGE, lw=2, label="Curva bolha (T–x)")
        ax2.plot(y, T, color=BLUE,   lw=2, ls="--", label="Curva orvalho (T–y)")
        ax2.fill_betweenx(T, x, y, alpha=0.08, color=GREEN)
        if az:
            ax2.axhline(az[2], color="red", lw=0.8, ls=":", alpha=0.7)
            ax2.axvline(az[0], color="red", lw=0.8, ls=":", alpha=0.7)
        ax2.legend(fontsize=8, facecolor=AXBG, labelcolor=TXT, edgecolor=GRID)

        # γ ou Psat
        ax3 = self.fig.add_subplot(gs[1,0])
        if nrtl:
            style_ax(ax3, "Coeficientes de atividade γ vs x₁", f"x₁  [{c1}]", "γ")
            ax3.plot(x, g1, color=ORANGE, lw=2, label=f"γ₁  {c1}")
            ax3.plot(x, g2, color=BLUE,   lw=2, label=f"γ₂  {c2}")
            ax3.axhline(1.0, color=GRAY, lw=1, ls=":", label="γ=1 (ideal)")
        else:
            style_ax(ax3, "Pressões de vapor vs Temperatura", "T (°C)", "Pˢᵃᵗ (bar)")
            ax3.plot(T, P1, color=ORANGE, lw=2, label=f"P₁ˢᵃᵗ  {c1}")
            ax3.plot(T, P2, color=BLUE,   lw=2, label=f"P₂ˢᵃᵗ  {c2}")
            ax3.axhline(P_atm, color=GRAY, lw=1, ls=":", label=f"P={P_atm:.3f} bar")
        ax3.legend(fontsize=8, facecolor=AXBG, labelcolor=TXT, edgecolor=GRID)

        # α₁₂
        ax4 = self.fig.add_subplot(gs[1,1])
        style_ax(ax4, "Volatilidade relativa α₁₂ vs x₁", f"x₁  [{c1}]", "α₁₂")
        ax4.plot(x, al, color=GREEN, lw=2)
        ax4.axhline(np.nanmean(al), color=GRAY, lw=1, ls="--",
                    label=f"α médio = {np.nanmean(al):.3f}")
        ax4.axhline(1.0, color=ORANGE, lw=0.8, ls=":", label="α=1")
        ax4.scatter(x[1:-1], al[1:-1], color=GREEN, s=18, zorder=5)
        if az:
            ax4.axvline(az[0], color="red", lw=0.8, ls=":", alpha=0.7,
                        label=f"Azeótropo x≈{az[0]:.3f}")
        ax4.legend(fontsize=8, facecolor=AXBG, labelcolor=TXT, edgecolor=GRID)

        self.fig.suptitle(
            f"Equilíbrio VLE — {c1} / {c2}   |   P={P_atm:.2f} atm   |   {lm} + Antoine (DDB/NIST)",
            color=TXT, fontsize=11, fontweight="bold", y=0.98)
        self.canvas.draw()

    # ── Tabela ────────────────────────────────
    def _preencher_tabela(self, df):
        for r in self.tree.get_children(): self.tree.delete(r)
        for _, r in df.iterrows():
            self.tree.insert("","end", values=(
                f"{r['x₁']:.4f}", f"{r['y₁']:.4f}",
                f"{r['T (K)']:.2f}", f"{r['T (°C)']:.2f}",
                f"{r['P₁ˢᵃᵗ (bar)']:.5f}", f"{r['P₂ˢᵃᵗ (bar)']:.5f}",
                f"{r['γ₁']:.4f}", f"{r['γ₂']:.4f}", f"{r['α₁₂ efetivo']:.4f}",
            ))

    # ── Exportação ────────────────────────────
    def _export_csv_full(self):
        if self.df_result is None:
            return messagebox.showwarning("Aviso","Calcule primeiro.")
        p = filedialog.asksaveasfilename(defaultextension=".csv",
                                          filetypes=[("CSV","*.csv")])
        if p:
            self.df_result.to_csv(p, index=False, sep=";", decimal=",")
            messagebox.showinfo("Exportado", f"Salvo em:\n{p}")

    def _export_csv_xy(self):
        if self.df_result is None:
            return messagebox.showwarning("Aviso","Calcule primeiro.")
        p = filedialog.asksaveasfilename(defaultextension=".csv",
                                          filetypes=[("CSV","*.csv")])
        if p:
            df_xy = self.df_result[["x₁","y₁"]].rename(columns={"x₁":"x","y₁":"y"})
            df_xy.to_csv(p, index=False, sep=",", decimal=".")
            messagebox.showinfo("Exportado", f"Salvo em:\n{p}")

    def _save_fig(self):
        p = filedialog.asksaveasfilename(defaultextension=".png",
                                          filetypes=[("PNG","*.png"),("PDF","*.pdf")])
        if p:
            self.fig.savefig(p, dpi=150, bbox_inches="tight", facecolor="#0E1117")
            messagebox.showinfo("Salvo", f"Salvo em:\n{p}")

    def _recarregar_banco(self):
        global COMPONENTS, NRTL_RAW, COMP_NAMES
        COMPONENTS  = carregar_componentes()
        NRTL_RAW    = carregar_nrtl()
        COMP_NAMES  = sorted(COMPONENTS.keys())
        self.ac1._names = COMP_NAMES
        self.ac2._names = COMP_NAMES
        messagebox.showinfo("Banco recarregado",
            f"✅ {len(COMPONENTS)} componentes e {len(NRTL_RAW)} pares NRTL carregados.")


if __name__ == "__main__":
    app = VLEApp()
    app.mainloop()
