"""
VLE Desktop — Equilíbrio Líquido-Vapor Bicomponente
Interface CustomTkinter + Matplotlib
Equivalente ao app Streamlit vle_app.py
"""

import math
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

# ─────────────────────────────────────────────
# Tema
# ─────────────────────────────────────────────
ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

# ─────────────────────────────────────────────
# Banco de componentes
# ─────────────────────────────────────────────
COMPONENT_GROUPS = {
    "── Aromáticos ──": None,
    "Benzeno":              {"A": 5.40768, "B": 1322.882, "C": -53.015,  "Tb_C": 80.1,   "MW": 78.11},
    "Tolueno":              {"A": 5.46600, "B": 1576.079, "C": -47.814,  "Tb_C": 110.6,  "MW": 92.14},
    "o-Xileno":             {"A": 5.51506, "B": 1736.072, "C": -46.900,  "Tb_C": 144.4,  "MW": 106.17},
    "m-Xileno":             {"A": 5.49750, "B": 1698.673, "C": -48.833,  "Tb_C": 139.1,  "MW": 106.17},
    "p-Xileno":             {"A": 5.49732, "B": 1691.879, "C": -49.235,  "Tb_C": 138.4,  "MW": 106.17},
    "Etilbenzeno":          {"A": 5.50675, "B": 1709.679, "C": -47.747,  "Tb_C": 136.2,  "MW": 106.17},
    "Estireno":             {"A": 5.56640, "B": 1861.894, "C": -44.500,  "Tb_C": 145.2,  "MW": 104.15},
    "── Alcanos ──": None,
    "n-Pentano":            {"A": 5.27087, "B": 1064.631, "C": -41.853,  "Tb_C": 36.1,   "MW": 72.15},
    "n-Hexano":             {"A": 5.26368, "B": 1202.948, "C": -52.636,  "Tb_C": 68.7,   "MW": 86.18},
    "n-Heptano":            {"A": 5.27786, "B": 1323.021, "C": -55.316,  "Tb_C": 98.4,   "MW": 100.20},
    "n-Octano":             {"A": 5.32054, "B": 1461.632, "C": -56.986,  "Tb_C": 125.7,  "MW": 114.23},
    "n-Nonano":             {"A": 5.33939, "B": 1575.415, "C": -60.410,  "Tb_C": 150.8,  "MW": 128.26},
    "Ciclohexano":          {"A": 5.26060, "B": 1295.030, "C": -58.100,  "Tb_C": 80.7,   "MW": 84.16},
    "Metilciclohexano":     {"A": 5.30550, "B": 1438.677, "C": -55.790,  "Tb_C": 100.9,  "MW": 98.19},
    "── Álcoois ──": None,
    "Metanol":              {"A": 5.31301, "B": 1676.569, "C": -21.728,  "Tb_C": 64.7,   "MW": 32.04},
    "Etanol":               {"A": 5.33675, "B": 1648.702, "C": -42.232,  "Tb_C": 78.4,   "MW": 46.07},
    "n-Propanol":           {"A": 5.37350, "B": 1788.020, "C": -35.940,  "Tb_C": 97.2,   "MW": 60.10},
    "i-Propanol":           {"A": 5.24268, "B": 1580.919, "C": -50.953,  "Tb_C": 82.4,   "MW": 60.10},
    "n-Butanol":            {"A": 5.36558, "B": 1891.523, "C": -36.055,  "Tb_C": 117.7,  "MW": 74.12},
    "── Cetonas e ésteres ──": None,
    "Acetona":              {"A": 5.31957, "B": 1490.864, "C": -35.930,  "Tb_C": 56.1,   "MW": 58.08},
    "MEK (butanona)":       {"A": 5.31424, "B": 1596.673, "C": -40.476,  "Tb_C": 79.6,   "MW": 72.11},
    "Acetato de etila":     {"A": 5.30680, "B": 1514.714, "C": -34.846,  "Tb_C": 77.1,   "MW": 88.11},
    "Acetato de n-butila":  {"A": 5.35647, "B": 1694.105, "C": -48.317,  "Tb_C": 126.1,  "MW": 116.16},
    "── Outros ──": None,
    "Água":                 {"A": 5.40221, "B": 1838.675, "C": -31.737,  "Tb_C": 100.0,  "MW": 18.02},
    "Clorofórmio":          {"A": 5.23628, "B": 1431.763, "C": -30.617,  "Tb_C": 61.2,   "MW": 119.38},
    "Diclorometano":        {"A": 5.20889, "B": 1325.938, "C": -24.064,  "Tb_C": 39.8,   "MW": 84.93},
    "Acetonitrila":         {"A": 5.28706, "B": 1492.380, "C": -32.613,  "Tb_C": 81.6,   "MW": 41.05},
    "Ácido acético":        {"A": 5.68206, "B": 1642.540, "C": -39.764,  "Tb_C": 117.9,  "MW": 60.05},
    "Dioxano":              {"A": 5.37096, "B": 1629.810, "C": -39.630,  "Tb_C": 101.3,  "MW": 88.11},
    "── Personalizado ──": None,
    "Personalizado":        {"A": None, "B": None, "C": None, "Tb_C": None, "MW": None},
}

COMPONENTS  = {k: v for k, v in COMPONENT_GROUPS.items() if v is not None}
COMP_OPTIONS = [k for k in COMPONENT_GROUPS.keys() if COMPONENT_GROUPS[k] is not None]

NRTL_PARAMS = {
    ("Etanol",  "Água"):                {"A12": 3458.3, "A21": -53.6,   "alpha": 0.30},
    ("Metanol", "Água"):                {"A12": 2726.4, "A21":  937.3,  "alpha": 0.30},
    ("Acetona", "Água"):                {"A12": 3768.3, "A21": 1116.0,  "alpha": 0.47},
    ("Acetona", "Clorofórmio"):         {"A12": -1704.0,"A21": -2051.0, "alpha": 0.30},
    ("Acetona", "Metanol"):             {"A12":  980.6, "A21":  481.2,  "alpha": 0.30},
    ("Etanol",  "Tolueno"):             {"A12": 4108.0, "A21": 1021.0,  "alpha": 0.47},
    ("Benzeno", "Etanol"):              {"A12": 3040.0, "A21": 1396.0,  "alpha": 0.47},
    ("Metanol", "Acetato de etila"):    {"A12": 1556.0, "A21": 1733.0,  "alpha": 0.47},
}

def get_nrtl_params(name1, name2):
    key, keyT = (name1, name2), (name2, name1)
    if key in NRTL_PARAMS:
        p = NRTL_PARAMS[key]
        return p["A12"], p["A21"], p["alpha"], False
    elif keyT in NRTL_PARAMS:
        p = NRTL_PARAMS[keyT]
        return p["A21"], p["A12"], p["alpha"], True
    return None, None, None, False

# ─────────────────────────────────────────────
# Termodinâmica
# ─────────────────────────────────────────────
def pvap(A, B, C, T_K):
    return 10.0 ** (A - B / (T_K + C))

def nrtl_gamma(x1, T_K, A12, A21, alpha):
    R = 8.314
    x2 = 1.0 - x1
    tau12 = A12 / (R * T_K)
    tau21 = A21 / (R * T_K)
    G12   = math.exp(-alpha * tau12)
    G21   = math.exp(-alpha * tau21)
    denom1 = x1 + x2 * G21
    denom2 = x2 + x1 * G12
    eps = 1e-12
    if abs(denom1) < eps: denom1 = eps
    if abs(denom2) < eps: denom2 = eps
    lng1 = x2**2 * (tau21 * (G21 / denom1)**2 + tau12 * G12 / denom2**2)
    lng2 = x1**2 * (tau12 * (G12 / denom2)**2 + tau21 * G21 / denom1**2)
    return math.exp(lng1), math.exp(lng2)

def bubble_T(x1, A1, B1, C1, A2, B2, C2, P_bar, T_init=None,
             modelo="Raoult", A12=0.0, A21=0.0, alpha_nrtl=0.30):
    if T_init is None:
        Tb1 = B1 / (A1 - math.log10(P_bar)) - C1
        Tb2 = B2 / (A2 - math.log10(P_bar)) - C2
        T   = x1 * Tb1 + (1.0 - x1) * Tb2
    else:
        T = T_init
    for _ in range(200):
        pb1   = pvap(A1, B1, C1, T)
        pb2   = pvap(A2, B2, C2, T)
        g1, g2 = (nrtl_gamma(x1, T, A12, A21, alpha_nrtl) if modelo == "NRTL" else (1.0, 1.0))
        Pcalc = x1 * g1 * pb1 + (1.0 - x1) * g2 * pb2
        err   = Pcalc - P_bar
        if abs(err) < 1e-7: break
        dT   = 0.01
        pb1d = pvap(A1, B1, C1, T + dT)
        pb2d = pvap(A2, B2, C2, T + dT)
        g1d, g2d = (nrtl_gamma(x1, T + dT, A12, A21, alpha_nrtl) if modelo == "NRTL" else (1.0, 1.0))
        dPdT = (x1 * g1d * pb1d + (1.0 - x1) * g2d * pb2d - Pcalc) / dT
        if abs(dPdT) < 1e-12: break
        T -= err / dPdT
    return T

def calc_vle(A1, B1, C1, A2, B2, C2, P_bar, n_points,
             modelo="Raoult", A12=0.0, A21=0.0, alpha_nrtl=0.30):
    x_arr  = np.linspace(0, 1, n_points)
    rows   = []
    T_prev = None
    for x1 in x_arr:
        T = bubble_T(x1, A1, B1, C1, A2, B2, C2, P_bar,
                     T_init=T_prev, modelo=modelo, A12=A12, A21=A21, alpha_nrtl=alpha_nrtl)
        T_prev = T
        P1s = pvap(A1, B1, C1, T)
        P2s = pvap(A2, B2, C2, T)
        g1, g2 = (nrtl_gamma(x1, T, A12, A21, alpha_nrtl) if modelo == "NRTL" else (1.0, 1.0))
        y1    = max(0.0, min(1.0, x1 * g1 * P1s / P_bar))
        alpha = (g1 * P1s) / (g2 * P2s) if P2s > 0 and g2 > 0 else np.nan
        rows.append({
            "x₁": round(float(x1), 6), "y₁": round(float(y1), 6),
            "T (K)": round(float(T), 4), "T (°C)": round(float(T - 273.15), 4),
            "P₁ˢᵃᵗ (bar)": round(float(P1s), 6), "P₂ˢᵃᵗ (bar)": round(float(P2s), 6),
            "γ₁": round(float(g1), 6), "γ₂": round(float(g2), 6),
            "α₁₂ efetivo": round(float(alpha), 4),
        })
    return pd.DataFrame(rows)

def detectar_azeotropo(df):
    diff = df["y₁"].values - df["x₁"].values
    for i in range(len(diff) - 1):
        if diff[i] * diff[i+1] < 0:
            x0, x1_ = df["x₁"].iloc[i], df["x₁"].iloc[i+1]
            d0, d1  = diff[i], diff[i+1]
            frac    = -d0 / (d1 - d0)
            x_az    = x0 + frac * (x1_ - x0)
            T_az    = df["T (°C)"].iloc[i] + frac * (df["T (°C)"].iloc[i+1] - df["T (°C)"].iloc[i])
            return x_az, x_az, T_az
    return None

# ─────────────────────────────────────────────
# Cores dos gráficos
# ─────────────────────────────────────────────
AXBG   = "#1A1D23"
GRID   = "#2A2D35"
TXT    = "#E0DDD6"
BLUE   = "#4A9EDF"
ORANGE = "#E87040"
GREEN  = "#5BB87A"
GRAY   = "#888880"

def style_ax(ax, title, xlabel, ylabel):
    ax.set_facecolor(AXBG)
    ax.tick_params(colors=TXT, labelsize=9)
    ax.xaxis.label.set_color(TXT)
    ax.yaxis.label.set_color(TXT)
    ax.title.set_color(TXT)
    ax.set_title(title, fontsize=10, fontweight="bold", pad=6)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.grid(color=GRID, linewidth=0.6, linestyle="--")
    for spine in ax.spines.values():
        spine.set_edgecolor(GRID)

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

        self.df_result   = None
        self.df_raoult   = None
        self._custom1    = {}
        self._custom2    = {}

        self._build_ui()
        self._on_comp1_change(self.comp1_var.get())
        self._on_comp2_change(self.comp2_var.get())
        self._on_modelo_change()
        self.after(100, self._calcular)

    # ── Layout principal ──────────────────────
    def _build_ui(self):
        # Coluna esquerda: painel de controle
        self.sidebar = ctk.CTkScrollableFrame(self, width=310, fg_color="#13161D",
                                               corner_radius=0)
        self.sidebar.pack(side="left", fill="y", padx=0, pady=0)

        # Área direita: gráficos + tabela
        self.main_area = ctk.CTkFrame(self, fg_color="#0E1117", corner_radius=0)
        self.main_area.pack(side="left", fill="both", expand=True)

        self._build_sidebar()
        self._build_main_area()

    def _section(self, parent, text):
        ctk.CTkLabel(parent, text=text, font=ctk.CTkFont(size=12, weight="bold"),
                     text_color="#4A9EDF").pack(anchor="w", padx=12, pady=(14, 2))
        ctk.CTkFrame(parent, height=1, fg_color="#2A2D35").pack(fill="x", padx=8, pady=(0, 8))

    def _label(self, parent, text):
        ctk.CTkLabel(parent, text=text, font=ctk.CTkFont(size=11),
                     text_color="#AAAAAA").pack(anchor="w", padx=14, pady=(4, 0))

    # ── Sidebar ───────────────────────────────
    def _build_sidebar(self):
        sb = self.sidebar

        # Título
        ctk.CTkLabel(sb, text="⚗️  VLE  Bicomponente",
                     font=ctk.CTkFont(size=15, weight="bold"),
                     text_color="#E0DDD6").pack(pady=(16, 4), padx=12, anchor="w")
        ctk.CTkLabel(sb, text="Raoult · NRTL · Antoine (NIST)",
                     font=ctk.CTkFont(size=10), text_color="#666666").pack(anchor="w", padx=14)

        # ── Modelo
        self._section(sb, "Modelo termodinâmico")
        self.modelo_var = ctk.StringVar(value="Raoult (ideal)")
        for opt in ["Raoult (ideal)", "NRTL (não ideal)"]:
            ctk.CTkRadioButton(sb, text=opt, variable=self.modelo_var, value=opt,
                               command=self._on_modelo_change,
                               font=ctk.CTkFont(size=11)).pack(anchor="w", padx=16, pady=2)

        # ── Comp 1
        self._section(sb, "Light Key (mais volátil)")
        self.comp1_var = ctk.StringVar(value="Benzeno")
        self.comp1_menu = ctk.CTkOptionMenu(sb, variable=self.comp1_var,
                                             values=COMP_OPTIONS, width=270,
                                             command=self._on_comp1_change)
        self.comp1_menu.pack(padx=12, pady=4)
        self.c1_info_label = ctk.CTkLabel(sb, text="", font=ctk.CTkFont(size=10),
                                           text_color="#666666")
        self.c1_info_label.pack(anchor="w", padx=14)

        # Campos custom comp1
        self.c1_custom_frame = ctk.CTkFrame(sb, fg_color="transparent")
        self._label(self.c1_custom_frame, "Nome")
        self.c1_nome = ctk.CTkEntry(self.c1_custom_frame, placeholder_text="Comp A", width=270)
        self.c1_nome.pack(padx=12, pady=2)
        for lbl, attr in [("A₁", "c1_A"), ("B₁", "c1_B"), ("C₁", "c1_C")]:
            self._label(self.c1_custom_frame, lbl)
            entry = ctk.CTkEntry(self.c1_custom_frame, width=270)
            entry.pack(padx=12, pady=2)
            setattr(self, attr, entry)

        # ── Comp 2
        self._section(sb, "Heavy Key (menos volátil)")
        self.comp2_var = ctk.StringVar(value="Tolueno")
        self.comp2_menu = ctk.CTkOptionMenu(sb, variable=self.comp2_var,
                                             values=COMP_OPTIONS, width=270,
                                             command=self._on_comp2_change)
        self.comp2_menu.pack(padx=12, pady=4)
        self.c2_info_label = ctk.CTkLabel(sb, text="", font=ctk.CTkFont(size=10),
                                           text_color="#666666")
        self.c2_info_label.pack(anchor="w", padx=14)

        self.c2_custom_frame = ctk.CTkFrame(sb, fg_color="transparent")
        self._label(self.c2_custom_frame, "Nome")
        self.c2_nome = ctk.CTkEntry(self.c2_custom_frame, placeholder_text="Comp B", width=270)
        self.c2_nome.pack(padx=12, pady=2)
        for lbl, attr in [("A₂", "c2_A"), ("B₂", "c2_B"), ("C₂", "c2_C")]:
            self._label(self.c2_custom_frame, lbl)
            entry = ctk.CTkEntry(self.c2_custom_frame, width=270)
            entry.pack(padx=12, pady=2)
            setattr(self, attr, entry)

        # ── NRTL params
        self.nrtl_frame = ctk.CTkFrame(sb, fg_color="transparent")
        self._section(self.nrtl_frame, "Parâmetros NRTL")
        self.nrtl_status = ctk.CTkLabel(self.nrtl_frame, text="", font=ctk.CTkFont(size=10),
                                         text_color="#5BB87A", wraplength=270)
        self.nrtl_status.pack(anchor="w", padx=14, pady=2)
        for lbl, attr, default in [("A₁₂ (J/mol)", "nrtl_A12", "1000.0"),
                                    ("A₂₁ (J/mol)", "nrtl_A21", "1000.0"),
                                    ("α (0.10 – 0.60)", "nrtl_alpha", "0.30")]:
            self._label(self.nrtl_frame, lbl)
            entry = ctk.CTkEntry(self.nrtl_frame, width=270)
            entry.insert(0, default)
            entry.pack(padx=12, pady=2)
            setattr(self, attr, entry)

        # ── Condições
        self._section(sb, "Condições de operação")
        self._label(sb, "Pressão (atm)")
        self.P_entry = ctk.CTkEntry(sb, width=270)
        self.P_entry.insert(0, "1.0")
        self.P_entry.pack(padx=12, pady=2)

        self._label(sb, "Número de pontos")
        self.n_slider_var = ctk.IntVar(value=21)
        self.n_slider = ctk.CTkSlider(sb, from_=11, to=101, number_of_steps=18,
                                       variable=self.n_slider_var, width=250,
                                       command=lambda v: self.n_label.configure(
                                           text=f"{int(float(v))} pontos"))
        self.n_slider.pack(padx=16, pady=4)
        self.n_label = ctk.CTkLabel(sb, text="21 pontos", font=ctk.CTkFont(size=10),
                                     text_color="#AAAAAA")
        self.n_label.pack(anchor="w", padx=14)

        # ── Botão calcular
        ctk.CTkButton(sb, text="🔄  Calcular", height=40, font=ctk.CTkFont(size=13, weight="bold"),
                       fg_color="#4A9EDF", hover_color="#3A8ECF",
                       command=self._calcular).pack(padx=12, pady=(16, 6), fill="x")

        # ── Botões exportar
        ctk.CTkButton(sb, text="⬇  Exportar CSV completo", height=34,
                       font=ctk.CTkFont(size=11), fg_color="#2A2D35", hover_color="#3A3D45",
                       command=self._export_csv_full).pack(padx=12, pady=3, fill="x")
        ctk.CTkButton(sb, text="⬇  Exportar x,y (McCabe-Thiele)", height=34,
                       font=ctk.CTkFont(size=11), fg_color="#2A2D35", hover_color="#3A3D45",
                       command=self._export_csv_xy).pack(padx=12, pady=3, fill="x")
        ctk.CTkButton(sb, text="🖼  Salvar gráficos (PNG)", height=34,
                       font=ctk.CTkFont(size=11), fg_color="#2A2D35", hover_color="#3A3D45",
                       command=self._save_fig).pack(padx=12, pady=(3, 16), fill="x")

    # ── Área principal ────────────────────────
    def _build_main_area(self):
        ma = self.main_area

        # Métricas no topo
        self.metrics_frame = ctk.CTkFrame(ma, height=60, fg_color="#13161D", corner_radius=8)
        self.metrics_frame.pack(fill="x", padx=12, pady=(10, 0))
        self.metrics_frame.pack_propagate(False)

        self.metric_labels = {}
        for i, key in enumerate(["P operação", "T_eb comp1", "T_eb comp2", "α₁₂ médio", "Azeótropo"]):
            f = ctk.CTkFrame(self.metrics_frame, fg_color="#1A1D23", corner_radius=6)
            f.pack(side="left", expand=True, fill="both", padx=6, pady=8)
            ctk.CTkLabel(f, text=key, font=ctk.CTkFont(size=9), text_color="#666666").pack(pady=(4, 0))
            lbl = ctk.CTkLabel(f, text="—", font=ctk.CTkFont(size=12, weight="bold"),
                                text_color="#E0DDD6")
            lbl.pack()
            self.metric_labels[key] = lbl

        # Abas: Gráficos | Tabela
        self.tabview = ctk.CTkTabview(ma, fg_color="#0E1117", segmented_button_fg_color="#13161D",
                                       segmented_button_selected_color="#4A9EDF",
                                       segmented_button_selected_hover_color="#3A8ECF",
                                       segmented_button_unselected_color="#13161D")
        self.tabview.pack(fill="both", expand=True, padx=12, pady=8)
        self.tab_graf = self.tabview.add("📈  Gráficos")
        self.tab_tab  = self.tabview.add("📋  Tabela")

        # Canvas matplotlib
        self.fig = plt.figure(figsize=(13, 8), facecolor="#0E1117")
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.tab_graf)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        # Tabela (Treeview via tkinter)
        self._build_table()

    def _build_table(self):
        frame = ctk.CTkFrame(self.tab_tab, fg_color="#0E1117")
        frame.pack(fill="both", expand=True, padx=8, pady=8)

        cols = ["x₁", "y₁", "T (K)", "T (°C)", "P₁ˢᵃᵗ (bar)", "P₂ˢᵃᵗ (bar)", "γ₁", "γ₂", "α₁₂ efetivo"]
        import tkinter.ttk as ttk

        style = ttk.Style()
        style.theme_use("default")
        style.configure("Treeview", background="#1A1D23", foreground="#E0DDD6",
                         fieldbackground="#1A1D23", rowheight=24, font=("Consolas", 10))
        style.configure("Treeview.Heading", background="#2A2D35", foreground="#4A9EDF",
                         font=("Consolas", 10, "bold"), relief="flat")
        style.map("Treeview", background=[("selected", "#4A9EDF")])

        self.tree = ttk.Treeview(frame, columns=cols, show="headings", height=30)
        for col in cols:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=110, anchor="center")

        vsb = ttk.Scrollbar(frame, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)

    # ── Callbacks ─────────────────────────────
    def _on_comp1_change(self, name):
        d = COMPONENTS.get(name, {})
        if name == "Personalizado":
            self.c1_custom_frame.pack(padx=0, pady=0, fill="x")
            self.c1_info_label.configure(text="")
        else:
            self.c1_custom_frame.pack_forget()
            if d:
                self.c1_info_label.configure(
                    text=f"A={d['A']}  B={d['B']}  C={d['C']}  |  T_eb={d['Tb_C']} °C")

    def _on_comp2_change(self, name):
        d = COMPONENTS.get(name, {})
        if name == "Personalizado":
            self.c2_custom_frame.pack(padx=0, pady=0, fill="x")
            self.c2_info_label.configure(text="")
        else:
            self.c2_custom_frame.pack_forget()
            if d:
                self.c2_info_label.configure(
                    text=f"A={d['A']}  B={d['B']}  C={d['C']}  |  T_eb={d['Tb_C']} °C")

    def _on_modelo_change(self):
        if self.modelo_var.get() == "NRTL (não ideal)":
            self.nrtl_frame.pack(fill="x")
            self._atualizar_nrtl_status()
        else:
            self.nrtl_frame.pack_forget()

    def _atualizar_nrtl_status(self):
        c1 = self.comp1_var.get() if self.comp1_var.get() != "Personalizado" else \
             (self.c1_nome.get() or "Comp A")
        c2 = self.comp2_var.get() if self.comp2_var.get() != "Personalizado" else \
             (self.c2_nome.get() or "Comp B")
        A12, A21, alpha, _ = get_nrtl_params(c1, c2)
        if A12 is not None:
            self.nrtl_status.configure(
                text=f"✅ Parâmetros encontrados: {c1} / {c2}", text_color="#5BB87A")
            self.nrtl_A12.delete(0, "end"); self.nrtl_A12.insert(0, str(A12))
            self.nrtl_A21.delete(0, "end"); self.nrtl_A21.insert(0, str(A21))
            self.nrtl_alpha.delete(0, "end"); self.nrtl_alpha.insert(0, str(alpha))
        else:
            self.nrtl_status.configure(
                text=f"⚠️ Par não encontrado — insira manualmente", text_color="#E87040")

    # ── Leitura dos parâmetros da UI ──────────
    def _get_params(self):
        c1_name = self.comp1_var.get()
        c2_name = self.comp2_var.get()

        if c1_name == "Personalizado":
            c1_label = self.c1_nome.get() or "Comp A"
            A1 = float(self.c1_A.get()); B1 = float(self.c1_B.get()); C1 = float(self.c1_C.get())
        else:
            d1 = COMPONENTS[c1_name]
            c1_label, A1, B1, C1 = c1_name, d1["A"], d1["B"], d1["C"]

        if c2_name == "Personalizado":
            c2_label = self.c2_nome.get() or "Comp B"
            A2 = float(self.c2_A.get()); B2 = float(self.c2_B.get()); C2 = float(self.c2_C.get())
        else:
            d2 = COMPONENTS[c2_name]
            c2_label, A2, B2, C2 = c2_name, d2["A"], d2["B"], d2["C"]

        P_atm = float(self.P_entry.get())
        n_pts = int(self.n_slider_var.get())
        usar_nrtl = self.modelo_var.get() == "NRTL (não ideal)"
        A12_val = float(self.nrtl_A12.get()) if usar_nrtl else 0.0
        A21_val = float(self.nrtl_A21.get()) if usar_nrtl else 0.0
        alpha_v = float(self.nrtl_alpha.get()) if usar_nrtl else 0.30

        return (c1_label, A1, B1, C1, c2_label, A2, B2, C2,
                P_atm, n_pts, usar_nrtl, A12_val, A21_val, alpha_v)

    # ── Cálculo principal ─────────────────────
    def _calcular(self):
        try:
            (c1_label, A1, B1, C1, c2_label, A2, B2, C2,
             P_atm, n_pts, usar_nrtl, A12_val, A21_val, alpha_v) = self._get_params()

            P_bar = P_atm

            # Verificar volatilidade
            Tb1 = bubble_T(1.0, A1, B1, C1, A2, B2, C2, P_bar)
            Tb2 = bubble_T(0.0, A1, B1, C1, A2, B2, C2, P_bar)
            swapped = False
            if Tb1 > Tb2:
                A1, B1, C1, A2, B2, C2 = A2, B2, C2, A1, B1, C1
                c1_label, c2_label     = c2_label, c1_label
                A12_val, A21_val       = A21_val, A12_val
                Tb1, Tb2               = Tb2, Tb1
                swapped = True

            modelo_str = "NRTL" if usar_nrtl else "Raoult"
            df = calc_vle(A1, B1, C1, A2, B2, C2, P_bar, n_pts,
                          modelo=modelo_str, A12=A12_val, A21=A21_val, alpha_nrtl=alpha_v)

            df_raoult = None
            if usar_nrtl:
                df_raoult = calc_vle(A1, B1, C1, A2, B2, C2, P_bar, n_pts, modelo="Raoult")

            az = detectar_azeotropo(df)
            self.df_result = df
            self.df_raoult = df_raoult
            self._last_params = (c1_label, c2_label, P_atm, usar_nrtl,
                                  A12_val, A21_val, alpha_v, az)

            # Métricas
            alpha_mean = df["α₁₂ efetivo"].mean()
            Teb1 = df["T (°C)"].iloc[-1]
            Teb2 = df["T (°C)"].iloc[0]
            self.metric_labels["P operação"].configure(text=f"{P_atm:.2f} atm")
            self.metric_labels["T_eb comp1"].configure(text=f"{Teb1:.1f} °C  ({c1_label})")
            self.metric_labels["T_eb comp2"].configure(text=f"{Teb2:.1f} °C  ({c2_label})")
            self.metric_labels["α₁₂ médio"].configure(text=f"{alpha_mean:.3f}")
            if az:
                self.metric_labels["Azeótropo"].configure(
                    text=f"x≈{az[0]:.3f}  T≈{az[2]:.1f}°C", text_color="#E87040")
            else:
                self.metric_labels["Azeótropo"].configure(text="Não detectado",
                                                           text_color="#5BB87A")

            self._plotar(df, df_raoult, c1_label, c2_label, P_atm, usar_nrtl, az)
            self._preencher_tabela(df)

            if swapped:
                messagebox.showinfo("Reordenação",
                    f"Componentes reordenados automaticamente:\n"
                    f"Light Key = {c1_label}\nHeavy Key = {c2_label}")
        except Exception as e:
            messagebox.showerror("Erro no cálculo", str(e))

    # ── Gráficos ──────────────────────────────
    def _plotar(self, df, df_raoult, c1, c2, P_atm, usar_nrtl, az):
        self.fig.clear()
        self.fig.patch.set_facecolor("#0E1117")
        gs = gridspec.GridSpec(2, 2, figure=self.fig, hspace=0.40, wspace=0.32,
                                left=0.07, right=0.97, top=0.93, bottom=0.07)

        x  = df["x₁"].values; y  = df["y₁"].values
        T  = df["T (°C)"].values
        P1 = df["P₁ˢᵃᵗ (bar)"].values; P2 = df["P₂ˢᵃᵗ (bar)"].values
        al = df["α₁₂ efetivo"].values
        g1 = df["γ₁"].values; g2 = df["γ₂"].values
        label_modelo = "NRTL" if usar_nrtl else "Raoult"

        # 1) y–x
        ax1 = self.fig.add_subplot(gs[0, 0])
        style_ax(ax1, f"Diagrama y–x  ({c1} / {c2})", f"x₁  [{c1}]", f"y₁  [{c1}]")
        if usar_nrtl and df_raoult is not None:
            ax1.plot(df_raoult["x₁"].values, df_raoult["y₁"].values,
                     color=BLUE, lw=1.5, ls="--", label="Raoult (ref.)", alpha=0.6)
        ax1.plot(x, y, color=ORANGE if usar_nrtl else BLUE, lw=2, label=label_modelo)
        ax1.plot([0, 1], [0, 1], color=GRAY, lw=1, ls="--", label="y = x")
        ax1.scatter(x[1:-1], y[1:-1], color=ORANGE if usar_nrtl else BLUE, s=18, zorder=5)
        if az:
            ax1.scatter([az[0]], [az[1]], color="red", s=80, zorder=10,
                        marker="*", label=f"Azeótropo x≈{az[0]:.3f}")
        ax1.set_xlim(0, 1); ax1.set_ylim(0, 1)
        ax1.legend(fontsize=8, facecolor=AXBG, labelcolor=TXT, edgecolor=GRID)

        # 2) T–x–y
        ax2 = self.fig.add_subplot(gs[0, 1])
        style_ax(ax2, f"Diagrama T–x–y  ({P_atm:.2f} atm)", f"Fração molar [{c1}]", "T (°C)")
        ax2.plot(x, T, color=ORANGE, lw=2, label="Curva bolha (T–x)")
        ax2.plot(y, T, color=BLUE,   lw=2, ls="--", label="Curva orvalho (T–y)")
        ax2.fill_betweenx(T, x, y, alpha=0.08, color=GREEN)
        if az:
            ax2.axhline(az[2], color="red", lw=0.8, ls=":", alpha=0.7)
            ax2.axvline(az[0], color="red", lw=0.8, ls=":", alpha=0.7)
        ax2.legend(fontsize=8, facecolor=AXBG, labelcolor=TXT, edgecolor=GRID)

        # 3) γ ou P_sat
        ax3 = self.fig.add_subplot(gs[1, 0])
        if usar_nrtl:
            style_ax(ax3, "Coeficientes de atividade γ vs x₁", f"x₁  [{c1}]", "γ")
            ax3.plot(x, g1, color=ORANGE, lw=2, label=f"γ₁  {c1}")
            ax3.plot(x, g2, color=BLUE,   lw=2, label=f"γ₂  {c2}")
            ax3.axhline(1.0, color=GRAY, lw=1, ls=":", label="γ = 1 (ideal)")
        else:
            style_ax(ax3, "Pressões de vapor vs Temperatura", "T (°C)", "Pˢᵃᵗ (bar)")
            ax3.plot(T, P1, color=ORANGE, lw=2, label=f"P₁ˢᵃᵗ  {c1}")
            ax3.plot(T, P2, color=BLUE,   lw=2, label=f"P₂ˢᵃᵗ  {c2}")
            ax3.axhline(float(self.P_entry.get()), color=GRAY, lw=1, ls=":",
                        label=f"P op = {float(self.P_entry.get()):.3f} bar")
        ax3.legend(fontsize=8, facecolor=AXBG, labelcolor=TXT, edgecolor=GRID)

        # 4) α₁₂
        ax4 = self.fig.add_subplot(gs[1, 1])
        style_ax(ax4, "Volatilidade relativa efetiva α₁₂ vs x₁", f"x₁  [{c1}]", "α₁₂")
        ax4.plot(x, al, color=GREEN, lw=2)
        ax4.axhline(np.nanmean(al), color=GRAY, lw=1, ls="--",
                    label=f"α médio = {np.nanmean(al):.3f}")
        ax4.axhline(1.0, color=ORANGE, lw=0.8, ls=":", label="α = 1")
        ax4.scatter(x[1:-1], al[1:-1], color=GREEN, s=18, zorder=5)
        if az:
            ax4.axvline(az[0], color="red", lw=0.8, ls=":", alpha=0.7,
                        label=f"Azeótropo x≈{az[0]:.3f}")
        ax4.legend(fontsize=8, facecolor=AXBG, labelcolor=TXT, edgecolor=GRID)

        self.fig.suptitle(
            f"Equilíbrio VLE — {c1} / {c2}   |   P = {P_atm:.2f} atm   |   {label_modelo} + Antoine (NIST)",
            color=TXT, fontsize=11, fontweight="bold", y=0.98)

        self.canvas.draw()

    # ── Tabela ────────────────────────────────
    def _preencher_tabela(self, df):
        for row in self.tree.get_children():
            self.tree.delete(row)
        for _, r in df.iterrows():
            self.tree.insert("", "end", values=(
                f"{r['x₁']:.4f}", f"{r['y₁']:.4f}",
                f"{r['T (K)']:.2f}", f"{r['T (°C)']:.2f}",
                f"{r['P₁ˢᵃᵗ (bar)']:.5f}", f"{r['P₂ˢᵃᵗ (bar)']:.5f}",
                f"{r['γ₁']:.4f}", f"{r['γ₂']:.4f}", f"{r['α₁₂ efetivo']:.4f}",
            ))

    # ── Exportação ────────────────────────────
    def _export_csv_full(self):
        if self.df_result is None:
            messagebox.showwarning("Aviso", "Calcule primeiro antes de exportar.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv",
                                            filetypes=[("CSV", "*.csv")],
                                            title="Salvar CSV completo")
        if path:
            self.df_result.to_csv(path, index=False, sep=";", decimal=",")
            messagebox.showinfo("Exportado", f"CSV salvo em:\n{path}")

    def _export_csv_xy(self):
        if self.df_result is None:
            messagebox.showwarning("Aviso", "Calcule primeiro antes de exportar.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv",
                                            filetypes=[("CSV", "*.csv")],
                                            title="Salvar x,y para McCabe-Thiele")
        if path:
            df_xy = self.df_result[["x₁", "y₁"]].rename(columns={"x₁": "x", "y₁": "y"})
            df_xy.to_csv(path, index=False, sep=",", decimal=".")
            messagebox.showinfo("Exportado", f"x,y salvo em:\n{path}")

    def _save_fig(self):
        path = filedialog.asksaveasfilename(defaultextension=".png",
                                            filetypes=[("PNG", "*.png"), ("PDF", "*.pdf")],
                                            title="Salvar gráficos")
        if path:
            self.fig.savefig(path, dpi=150, bbox_inches="tight",
                             facecolor="#0E1117")
            messagebox.showinfo("Salvo", f"Gráfico salvo em:\n{path}")


# ─────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────
if __name__ == "__main__":
    app = VLEApp()
    app.mainloop()
