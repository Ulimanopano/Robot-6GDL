# -*- coding: utf-8 -*-
"""
Robot DH Modificado | Editor de parámetros y cinemática directa.

Requiere Python >= 3.10 y SymPy:
    py -m pip install sympy
Ejecutar:
    py Robot_DH_Modificado.py

Convención MDH (Craig):
    T_(i-1,i) = Rx(alpha_(i-1)) Tx(a_(i-1)) Rz(theta_i) Tz(d_i)
    Ejemplo: theta2-pi/2, -pi/2, 30°, a1, L1, Lef.

Comprobación sin interfaz:
    py Robot_DH_Modificado.py --self-test
"""

from __future__ import annotations

import ast
import json
import queue
import re
import sys
import threading
from pathlib import Path

try:
    import sympy as sp
except ImportError:
    print("Falta SymPy. Instálalo con: py -m pip install sympy")
    raise SystemExit(1)


# ========================== CÁLCULO SIMBÓLICO ===============================

ALLOWED_FUNCTIONS = {
    "sin": sp.sin, "cos": sp.cos, "tan": sp.tan,
    "sqrt": sp.sqrt, "asin": sp.asin, "acos": sp.acos,
    "atan": sp.atan, "abs": sp.Abs,
}
CONSTANTS = {"pi": sp.pi, "E": sp.E, "tau": 2 * sp.pi}


def parse_expression(source: str) -> sp.Expr:
    """Convierte expresiones del usuario a SymPy, sin ejecutar código Python.

    Admite: +, -, *, /, ^, (), pi, θ, π, grados 30° y funciones básicas.
    No admite multiplicación implícita (escribir 2*a1, no 2a1).
    """
    source = str(source).strip()
    if not source:
        raise ValueError("Hay una celda vacía. Usa 0 si el parámetro es nulo.")
    if len(source) > 180:
        raise ValueError("La expresión es demasiado larga (máximo 180 caracteres).")
    source = (source.replace("π", "pi")
                    .replace("θ", "theta")
                    .replace("−", "-")
                    .replace("×", "*")
                    .replace("÷", "/")
                    .replace("^", "**"))
    # El sufijo de grados funciona en el campo de ángulo y en los valores.
    source = re.sub(
        r"(?<![\w.])(\d+(?:\.\d+)?)\s*(?:°|\bgrados\b|\bdeg\b)",
        r"(\1*pi/180)", source, flags=re.IGNORECASE
    )

    try:
        parsed = ast.parse(source, mode="eval")
    except SyntaxError as ex:
        raise ValueError(f"Expresión incorrecta: {source}") from ex

    def convert(node: ast.AST) -> sp.Expr:
        if isinstance(node, ast.Constant):
            val = node.value
            if isinstance(val, bool) or not isinstance(val, (int, float)):
                raise ValueError("Solo se admiten números y variables simbólicas.")
            return sp.Rational(str(val))
        if isinstance(node, ast.Name):
            if node.id in CONSTANTS:
                return CONSTANTS[node.id]
            if node.id in ALLOWED_FUNCTIONS:
                raise ValueError("A una función debes agregarle paréntesis, ej. sin(theta1).")
            if node.id.startswith("__"):
                raise ValueError("Nombre de variable no permitido.")
            return sp.Symbol(node.id, real=True)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
            sub = convert(node.operand)
            return -sub if isinstance(node.op, ast.USub) else sub
        if isinstance(node, ast.BinOp):
            left, right = convert(node.left), convert(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
            if isinstance(node.op, ast.Pow):
                if not right.is_number or abs(float(right.evalf())) > 8:
                    raise ValueError("Solo se admiten potencias constantes de -8 a 8.")
                return left ** right
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            name = node.func.id
            if name not in ALLOWED_FUNCTIONS or len(node.args) != 1 or node.keywords:
                raise ValueError("Función no admitida. Usa sin, cos, tan, sqrt, asin, acos, atan o abs.")
            return ALLOWED_FUNCTIONS[name](convert(node.args[0]))
        raise ValueError("Sintaxis no permitida. Usa números, variables, + - * / y paréntesis.")

    try:
        out = convert(parsed.body)
    except (OverflowError, ZeroDivisionError, TypeError) as ex:
        raise ValueError(f"No se pudo interpretar {source!r}: {ex}") from ex
    if out.has(sp.zoo, sp.nan, sp.oo, -sp.oo):
        raise ValueError("La expresión contiene una división entre cero o un infinito.")
    return out


def mdh_transform(alpha: sp.Expr, a: sp.Expr, d: sp.Expr, theta: sp.Expr) -> sp.Matrix:
    """DH modificada de Craig: Rx(alpha) Tx(a) Rz(theta) Tz(d)."""
    ca, sa = sp.cos(alpha), sp.sin(alpha)
    ct, st = sp.cos(theta), sp.sin(theta)
    return sp.Matrix([
        [ct, -st, 0, a],
        [st * ca, ct * ca, -sa, -d * sa],
        [st * sa, ct * sa, ca, d * ca],
        [0, 0, 0, 1],
    ])


def compute_kinematics(rows: list[dict[str, str]], values: dict[str, str] | None = None,
                       simplify: bool = False) -> tuple[list[sp.Matrix], list[sp.Matrix]]:
    """Devuelve ([T01,T12,...], [T01,T02,...,T0n]).

    Valores ausentes/vacíos => variables simbólicas; valores presentes => constantes.
    Los ángulos sin ° se interpretan en radianes.
    """
    if not rows:
        raise ValueError("La tabla DH debe contener al menos una fila.")
    parsed_rows = []
    for i, row in enumerate(rows, start=1):
        try:
            parsed_rows.append({key: parse_expression(row[key]) for key in ("alpha", "a", "d", "theta")})
        except (ValueError, KeyError) as ex:
            raise ValueError(f"Error en fila {i}: {ex}") from ex

    substitutions = {}
    for name, string_val in (values or {}).items():
        if not str(string_val).strip():
            continue
        val = parse_expression(string_val)
        if val.free_symbols:
            raise ValueError(f"El valor de {name} debe ser constante; contiene {val.free_symbols}.")
        substitutions[sp.Symbol(name, real=True)] = val

    individual = []
    cumulative = []
    total = sp.eye(4)
    for row in parsed_rows:
        params = [row[key].subs(substitutions) for key in ("alpha", "a", "d", "theta")]
        transform = mdh_transform(*params)
        total = total * transform
        if simplify:
            # Puede tardar con 6+ GDL; por eso se ejecuta en segundo plano.
            transform = transform.applyfunc(sp.trigsimp)
            total = total.applyfunc(sp.trigsimp)
        individual.append(transform)
        cumulative.append(total)
    return individual, cumulative


def natural_sort_key(value: str):
    """Agrupa theta1,...,theta10 de manera intuitiva."""
    groups = re.split(r"(\d+)", value.lower())
    chunks = tuple(int(chunk) if chunk.isdigit() else chunk for chunk in groups)
    return (0 if value.lower().startswith("theta") else 1, chunks)


def pretty_matrix(mat: sp.Matrix) -> str:
    return sp.pretty(mat, use_unicode=True, wrap_line=False)


def example_rows() -> list[dict[str, str]]:
    """Última fotografía enviada: 3 GDL y un eslabón herramienta fijo."""
    return [
        {"alpha": "0", "a": "0", "d": "L1", "theta": "theta1"},
        {"alpha": "-pi/2", "a": "a1", "d": "L2", "theta": "theta2-pi/2"},
        {"alpha": "-pi/2", "a": "0", "d": "-a2", "theta": "theta3-pi/2"},
        {"alpha": "0", "a": "Lef", "d": "0", "theta": "0"},
    ]


def self_test():
    """Verifica el ejemplo de clase y el comportamiento simbólico."""
    rows = example_rows()
    consecutive, accumulated = compute_kinematics(rows)
    assert len(consecutive) == len(accumulated) == 4
    assert all(M.shape == (4, 4) for M in consecutive + accumulated)
    theta = {f"theta{i}": "0" for i in range(1, 4)}
    _, evaluated = compute_kinematics(rows, theta)
    a1, a2, L1, L2, Lef = sp.symbols("a1 a2 L1 L2 Lef", real=True)
    expected = sp.Matrix([
        [0, 0, 1, a1 - a2],
        [1, 0, 0, L2 + Lef],
        [0, 1, 0, L1],
        [0, 0, 0, 1],
    ])
    assert evaluated[-1] == expected, f"Esperado:\n{expected}\nObtenido:\n{evaluated[-1]}"
    assert parse_expression("30°") == sp.pi / 6
    assert parse_expression("θ3-π/2") == sp.Symbol("theta3", real=True) - sp.pi / 2
    assert parse_expression("3*a1") == 3 * sp.Symbol("a1", real=True)
    assert compute_kinematics([{"alpha":"0","a":"0","d":"0","theta":"theta1"}])[0][0][0,0] == sp.cos(sp.Symbol("theta1", real=True))
    print("PRUEBAS CORRECTAS\n")
    print("Última tabla de tu clase, con theta1=theta2=theta3=0:\n")
    print(pretty_matrix(evaluated[-1]))


# ============================== INTERFAZ ====================================


def start_gui():
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox, ttk
    except ImportError:
        print("No se encontró Tkinter. Instala una distribución de Python que incluya Tk/Tcl.")
        raise SystemExit(1)

    class RobotDHApp:
        def __init__(self, root):
            self.root = root
            self.root.title("Robot DH Modificado  |  Cinemática directa simbólica")
            self.root.geometry("1220x850")
            self.root.minsize(930, 650)
            self._building = False
            self._debounce_id = None
            self._in_progress = False
            self.prev_gdl = 3
            self.prev_has_tool = True
            self.rows: list[dict[str, tk.StringVar]] = []
            self.variable_widgets: dict[str, tk.StringVar] = {}
            self.saved_values: dict[str, str] = {}
            self.last_matrices = None
            self.last_output_text = {}
            self.events = queue.Queue()
            self.gdl = tk.StringVar(value="3")
            self.has_tool = tk.BooleanVar(value=True)
            self.simplify = tk.BooleanVar(value=False)
            self.status = tk.StringVar(value="Listo. Edita las celdas y pulsa «Calcular matrices».")

            self._style()
            self._layout()
            self._render_table(example_rows())
            self._refresh_variables()
            self._set_tab_text("Final T0n", "Pulsa «Calcular matrices» para ver la transformación total.")

        def _style(self):
            style = ttk.Style(self.root)
            try:
                style.theme_use("clam")
            except tk.TclError:
                pass
            style.configure("TButton", font=("Segoe UI", 10), padding=(10, 7))
            style.configure("Primary.TButton", font=("Segoe UI", 10, "bold"), padding=(14, 8))
            style.configure("TLabel", font=("Segoe UI", 10))
            style.configure("Title.TLabel", font=("Segoe UI", 18, "bold"), foreground="#24324b")
            style.configure("Info.TLabel", font=("Segoe UI", 9), foreground="#57657a")
            style.configure("TLabelframe.Label", font=("Segoe UI", 10, "bold"))
            style.configure("TNotebook.Tab", font=("Segoe UI", 10), padding=(13, 7))

        def _layout(self):
            base = ttk.Frame(self.root, padding=14)
            base.pack(fill="both", expand=True)

            ttk.Label(base, text="Cinemática directa · Denavit–Hartenberg modificada", style="Title.TLabel").pack(anchor="w")
            ttk.Label(base,
                text="Igual a la tabla de tu pizarrón: α(i−1), a(i−1), d(i), θ(i)  |  Rx(α) · Tx(a) · Rz(θ) · Tz(d)",
                style="Info.TLabel").pack(anchor="w", pady=(2, 12))

            control = ttk.Frame(base)
            control.pack(fill="x", pady=(0, 10))
            ttk.Label(control, text="Articulaciones / GDL:").pack(side="left", padx=(0, 7))
            self.gdl_spin = ttk.Spinbox(control, from_=1, to=12, width=5, textvariable=self.gdl,
                                        command=self._change_structure)
            self.gdl_spin.pack(side="left")
            self.gdl_spin.bind("<Return>", lambda e: self._change_structure())
            self.gdl_spin.bind("<FocusOut>", lambda e: self._change_structure())
            ttk.Checkbutton(control, text="Añadir efector / herramienta fija", variable=self.has_tool,
                            command=self._change_structure).pack(side="left", padx=(18, 16))
            ttk.Checkbutton(control, text="Simplificar trigonometría (más lento)", variable=self.simplify).pack(side="left")
            ttk.Button(control, text="Ejemplo del pizarrón", command=self._load_example).pack(side="right")

            config = ttk.Panedwindow(base, orient="horizontal")
            config.pack(fill="both", expand=True)
            params_outer = ttk.LabelFrame(config, text="  Tabla de parámetros DH modificada  ")
            variables_outer = ttk.LabelFrame(config, text="  Valores de las variables  ")
            config.add(params_outer, weight=4)
            config.add(variables_outer, weight=2)

            self.param_canvas = tk.Canvas(params_outer, borderwidth=0, highlightthickness=0, height=225)
            p_scroll = ttk.Scrollbar(params_outer, orient="vertical", command=self.param_canvas.yview)
            self.param_canvas.configure(yscrollcommand=p_scroll.set)
            self.param_canvas.pack(side="left", fill="both", expand=True)
            p_scroll.pack(side="right", fill="y")
            self.param_body = ttk.Frame(self.param_canvas, padding=8)
            self.param_window = self.param_canvas.create_window((0, 0), window=self.param_body, anchor="nw")
            self.param_body.bind("<Configure>", self._update_param_scroll)
            self.param_canvas.bind("<Configure>", self._resize_param_width)

            variables_note = ttk.Label(variables_outer, text="Vacío = simbólica   |   45° = grados   |   pi/4 = radianes",
                                       style="Info.TLabel", wraplength=350)
            variables_note.pack(fill="x", padx=9, pady=(9, 5))
            self.var_canvas = tk.Canvas(variables_outer, borderwidth=0, highlightthickness=0, height=200)
            v_scroll = ttk.Scrollbar(variables_outer, orient="vertical", command=self.var_canvas.yview)
            self.var_canvas.configure(yscrollcommand=v_scroll.set)
            self.var_canvas.pack(side="left", fill="both", expand=True)
            v_scroll.pack(side="right", fill="y")
            self.var_body = ttk.Frame(self.var_canvas, padding=8)
            self.var_window = self.var_canvas.create_window((0, 0), window=self.var_body, anchor="nw")
            self.var_body.bind("<Configure>", lambda _: self.var_canvas.configure(scrollregion=self.var_canvas.bbox("all")))
            self.var_canvas.bind("<Configure>", lambda e: self.var_canvas.itemconfigure(self.var_window, width=e.width))

            example = ttk.Label(base,
                text="Ejemplos: θ = theta2-pi/2, d = -a2, a = Lef, α = -pi/2. Para multiplicar usa 2*a1. "
                     "Un ángulo numérico sin ° está en radianes.",
                style="Info.TLabel", wraplength=1150)
            example.pack(fill="x", pady=(7, 6))

            actions = ttk.Frame(base)
            actions.pack(fill="x", pady=(2, 9))
            self.calculate_button = ttk.Button(actions, text="CALCULAR MATRICES", style="Primary.TButton",
                                               command=self._calculate)
            self.calculate_button.pack(side="left", padx=(0, 8))
            ttk.Button(actions, text="Guardar proyecto", command=self._save_project).pack(side="left", padx=4)
            ttk.Button(actions, text="Abrir proyecto", command=self._open_project).pack(side="left", padx=4)
            ttk.Button(actions, text="Exportar TXT", command=self._export_txt).pack(side="right", padx=3)
            ttk.Button(actions, text="Exportar LaTeX", command=self._export_tex).pack(side="right", padx=3)
            ttk.Button(actions, text="Copiar pestaña", command=self._copy_tab).pack(side="right", padx=3)

            output = ttk.LabelFrame(base, text="  Resultados  ")
            output.pack(fill="both", expand=True)
            self.notebook = ttk.Notebook(output)
            self.notebook.pack(fill="both", expand=True, padx=6, pady=6)
            self.result_texts = {}
            for label in ("Final T0n", "Consecutivas", "Acumuladas"):
                frame = ttk.Frame(self.notebook)
                self.notebook.add(frame, text=label)
                result = tk.Text(frame, wrap="none", font=("Consolas", 12), bg="#fafbfe",
                                 foreground="#263547", undo=False, relief="flat", padx=12, pady=12)
                y = ttk.Scrollbar(frame, orient="vertical", command=result.yview)
                x = ttk.Scrollbar(frame, orient="horizontal", command=result.xview)
                result.configure(yscrollcommand=y.set, xscrollcommand=x.set)
                x.pack(side="bottom", fill="x")
                y.pack(side="right", fill="y")
                result.pack(side="left", fill="both", expand=True)
                self.result_texts[label] = result
            ttk.Label(base, textvariable=self.status, style="Info.TLabel").pack(fill="x", pady=(7, 0))

        def _update_param_scroll(self, _):
            self.param_canvas.configure(scrollregion=self.param_canvas.bbox("all"))

        def _resize_param_width(self, event):
            self.param_canvas.itemconfigure(self.param_window, width=event.width)

        def _get_gdl(self):
            try:
                v = int(self.gdl.get())
            except ValueError as ex:
                raise ValueError("Introduce un número entero de GDL entre 1 y 12.") from ex
            if not 1 <= v <= 12:
                raise ValueError("Los GDL deben estar entre 1 y 12.")
            return v

        def _render_table(self, source_rows: list[dict[str, str]]):
            self._building = True
            for child in self.param_body.winfo_children():
                child.destroy()
            headers = ["i", "α (i−1)", "a (i−1)", "d (i)", "θ (i)"]
            for j, h in enumerate(headers):
                ttk.Label(self.param_body, text=h, font=("Segoe UI", 10, "bold"),
                          foreground="#2e4c75").grid(row=0, column=j, padx=4, pady=(2, 11), sticky="w")
            self.rows = []
            joints = self._get_gdl()
            for i, source in enumerate(source_rows, start=1):
                is_tool = self.has_tool.get() and i == joints + 1
                text = f"{i}  · herramienta" if is_tool else str(i)
                ttk.Label(self.param_body, text=text).grid(row=i, column=0, padx=4, pady=6, sticky="w")
                record = {}
                for j, key in enumerate(("alpha", "a", "d", "theta"), start=1):
                    value = tk.StringVar(value=source.get(key, "0"))
                    entry = ttk.Entry(self.param_body, textvariable=value, width=17,
                                      font=("Consolas", 11))
                    entry.grid(row=i, column=j, sticky="ew", padx=4, pady=6)
                    value.trace_add("write", lambda *_: self._schedule_variable_sync())
                    record[key] = value
                self.rows.append(record)
            for j in range(1, 5):
                self.param_body.columnconfigure(j, weight=1)
            self._building = False
            self.param_body.update_idletasks()
            self._refresh_variables()

        def _snapshot_rows(self):
            return [{key: record[key].get() for key in ("alpha", "a", "d", "theta")}
                    for record in self.rows]

        def _new_row(self, i, tool=False):
            return {"alpha": "0", "a": "Lef" if tool else "0", "d": "0",
                    "theta": "0" if tool else f"theta{i}"}

        def _change_structure(self):
            if self._building:
                return
            try:
                new_gdl = self._get_gdl()
            except ValueError as ex:
                self.status.set(str(ex))
                return
            new_tool = bool(self.has_tool.get())
            if new_gdl == self.prev_gdl and new_tool == self.prev_has_tool:
                return
            previous = self._snapshot_rows()
            old_joint_rows = previous[:self.prev_gdl]
            old_tool_row = previous[self.prev_gdl] if self.prev_has_tool and len(previous) > self.prev_gdl else None
            new_rows = [old_joint_rows[i] if i < len(old_joint_rows) else self._new_row(i+1)
                        for i in range(new_gdl)]
            if new_tool:
                new_rows.append(old_tool_row or self._new_row(new_gdl+1, tool=True))
            self.prev_gdl, self.prev_has_tool = new_gdl, new_tool
            self._render_table(new_rows)
            self.status.set(f"Tabla preparada: {new_gdl} articulaciones y {len(new_rows)} transformaciones.")

        def _schedule_variable_sync(self):
            if self._building:
                return
            if self._debounce_id is not None:
                self.root.after_cancel(self._debounce_id)
            self._debounce_id = self.root.after(360, self._refresh_variables)

        def _snapshot_values(self):
            stored = dict(self.saved_values)
            for name, var in self.variable_widgets.items():
                stored[name] = var.get()
            return stored

        def _refresh_variables(self):
            if self._building:
                return
            self._debounce_id = None
            last_values = self._snapshot_values()
            symbols = set()
            errors = []
            for i, row in enumerate(self._snapshot_rows(), start=1):
                for field, source in row.items():
                    try:
                        symbols.update(str(sym) for sym in parse_expression(source).free_symbols)
                    except ValueError:
                        errors.append(f"Fila {i}, {field}")
            for child in self.var_body.winfo_children():
                child.destroy()
            self.variable_widgets = {}
            if not symbols:
                ttk.Label(self.var_body, text="No hay variables simbólicas.",
                          style="Info.TLabel").grid(row=0, column=0, sticky="w")
            else:
                for i, name in enumerate(sorted(symbols, key=natural_sort_key)):
                    ttk.Label(self.var_body, text=f"{name} =", font=("Consolas", 11)).grid(
                        row=i, column=0, sticky="w", padx=(1, 9), pady=5)
                    var = tk.StringVar(value=last_values.get(name, ""))
                    ttk.Entry(self.var_body, textvariable=var, width=18, font=("Consolas", 11)).grid(
                        row=i, column=1, sticky="ew", padx=(0, 6), pady=5)
                    self.variable_widgets[name] = var
                self.var_body.columnconfigure(1, weight=1)
            self.saved_values = last_values
            if errors:
                self.status.set("Revisa expresiones incompletas en: " + ", ".join(errors[:5]))
            else:
                self.status.set(f"{len(symbols)} variables reconocidas automáticamente. Deja en blanco las simbólicas.")

        def _load_example(self):
            self.gdl.set("3")
            self.has_tool.set(True)
            self.prev_gdl = 3
            self.prev_has_tool = True
            self.saved_values = {}
            self.variable_widgets = {}
            self._render_table(example_rows())
            self.status.set("Ejemplo de tu última fotografía: 3 GDL + herramienta fija.")

        def _calculate(self):
            if self._in_progress:
                return
            try:
                n = self._get_gdl()
                if n != self.prev_gdl or self.has_tool.get() != self.prev_has_tool:
                    self._change_structure()
                rows = self._snapshot_rows()
                values = self._snapshot_values()
                # Validación inmediata en el hilo principal (rápida).
                for idx, row in enumerate(rows, start=1):
                    for key, source in row.items():
                        try:
                            parse_expression(source)
                        except ValueError as ex:
                            raise ValueError(f"Fila {idx}, {key}: {ex}") from ex
                symbols = set().union(*(parse_expression(x).free_symbols for row in rows for x in row.values()))
                values = {str(s): values.get(str(s), "") for s in symbols}
                for sym, source in values.items():
                    if source.strip() and parse_expression(source).free_symbols:
                        raise ValueError(f"La variable {sym} necesita un valor numérico, no otra variable.")
            except ValueError as ex:
                messagebox.showerror("Error en los parámetros", str(ex))
                return
            self._in_progress = True
            self.calculate_button.state(["disabled"])
            self.status.set("Calculando transformaciones simbólicas... (6+ GDL pueden tardar)")
            worker = threading.Thread(
                target=self._worker, args=(rows, values, self.simplify.get()), daemon=True)
            worker.start()
            self.root.after(100, self._poll)

        def _worker(self, rows, values, simplify):
            try:
                individual, cumulative = compute_kinematics(rows, values, simplify)
                count = len(individual)
                consecutive_output = [
                    f"T({i-1},{i}) =\n{pretty_matrix(mat)}"
                    for i, mat in enumerate(individual, start=1)
                ]
                accumulated_output = [
                    f"T(0,{i}) =\n{pretty_matrix(mat)}"
                    for i, mat in enumerate(cumulative, start=1)
                ]
                result = {
                    "Final T0n": f"TRANSFORMACIÓN TOTAL T(0,{count})\n\n"
                                   + pretty_matrix(cumulative[-1])
                                   + "\n\nPOSICIÓN DEL EFECTOR (respecto al marco 0):\n"
                                   + pretty_matrix(cumulative[-1][:3, 3])
                                   + "\n\nORIENTACIÓN R(0,n):\n"
                                   + pretty_matrix(cumulative[-1][:3, :3]),
                    "Consecutivas": "MATRICES ENTRE SISTEMAS CONSECUTIVOS\n\n" + "\n\n".join(consecutive_output),
                    "Acumuladas": "MATRICES DESDE EL SISTEMA BASE 0\n\n" + "\n\n".join(accumulated_output),
                }
                self.events.put((True, result, individual, cumulative))
            except Exception as ex:
                self.events.put((False, str(ex), None, None))

        def _poll(self):
            try:
                ok, payload, individual, cumulative = self.events.get_nowait()
            except queue.Empty:
                self.root.after(100, self._poll)
                return
            self._in_progress = False
            self.calculate_button.state(["!disabled"])
            if not ok:
                self.status.set("Ocurrió un error en el cálculo.")
                messagebox.showerror("No se pudieron calcular las matrices", payload)
                return
            self.last_matrices = (individual, cumulative)
            self.last_output_text = payload
            for label, content in payload.items():
                self._set_tab_text(label, content)
            self.status.set(f"Completado: {len(individual)} consecutivas, "
                            f"{len(cumulative)} acumuladas y matriz final T(0,{len(cumulative)}).")

        def _set_tab_text(self, label, content):
            text_widget = self.result_texts[label]
            text_widget.configure(state="normal")
            text_widget.delete("1.0", "end")
            text_widget.insert("end", content)
            text_widget.configure(state="disabled")

        def _current_tab_name(self):
            return self.notebook.tab(self.notebook.select(), "text")

        def _copy_tab(self):
            name = self._current_tab_name()
            content = self.result_texts[name].get("1.0", "end-1c")
            self.root.clipboard_clear()
            self.root.clipboard_append(content)
            self.status.set(f"Contenido de «{name}» copiado al portapapeles.")

        def _save_project(self):
            try:
                self._get_gdl()
                data = {
                    "format": "RobotDHModificado-1",
                    "gdl": self.prev_gdl,
                    "herramienta_fija": self.prev_has_tool,
                    "simplificar": self.simplify.get(),
                    "rows": self._snapshot_rows(),
                    "values": self._snapshot_values(),
                }
                path = filedialog.asksaveasfilename(title="Guardar proyecto DH", defaultextension=".json",
                                                    filetypes=[("Proyecto DH", "*.json")])
                if path:
                    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
                    self.status.set(f"Proyecto guardado: {path}")
            except (OSError, ValueError) as ex:
                messagebox.showerror("No se pudo guardar", str(ex))

        def _open_project(self):
            path = filedialog.askopenfilename(title="Abrir proyecto DH", filetypes=[("Proyecto DH", "*.json")])
            if not path:
                return
            try:
                data = json.loads(Path(path).read_text(encoding="utf-8"))
                if data.get("format") != "RobotDHModificado-1":
                    raise ValueError("Archivo JSON incompatible con este programa.")
                n = int(data["gdl"])
                tool = bool(data["herramienta_fija"])
                if not 1 <= n <= 12 or len(data["rows"]) != n + int(tool):
                    raise ValueError("Cantidad de filas / GDL inconsistente.")
                self.gdl.set(str(n))
                self.has_tool.set(tool)
                self.prev_gdl, self.prev_has_tool = n, tool
                self.simplify.set(bool(data.get("simplificar", False)))
                self.saved_values = dict(data.get("values", {}))
                self.variable_widgets = {}
                self._render_table(data["rows"])
                self.status.set(f"Proyecto cargado: {path}")
            except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as ex:
                messagebox.showerror("No se pudo abrir el proyecto", str(ex))

        def _export_txt(self):
            if not self.last_output_text:
                messagebox.showinfo("Sin resultados", "Primero pulsa «Calcular matrices».")
                return
            path = filedialog.asksaveasfilename(title="Exportar resultados de matrices", defaultextension=".txt",
                                                filetypes=[("Texto", "*.txt")])
            if not path:
                return
            text = "CINEMÁTICA DIRECTA — DH MODIFICADA\n" + "=" * 54 + "\n\n"
            for tab in ("Final T0n", "Consecutivas", "Acumuladas"):
                text += f"\n{'-'*24} {tab} {'-'*24}\n\n{self.last_output_text[tab]}\n"
            try:
                Path(path).write_text(text, encoding="utf-8")
                self.status.set(f"Matrices exportadas a {path}")
            except OSError as ex:
                messagebox.showerror("No se pudo exportar", str(ex))

        def _export_tex(self):
            if not self.last_matrices:
                messagebox.showinfo("Sin resultados", "Primero pulsa «Calcular matrices».")
                return
            path = filedialog.asksaveasfilename(title="Exportar matrices en LaTeX", defaultextension=".tex",
                                                filetypes=[("LaTeX", "*.tex")])
            if not path:
                return
            individual, cumulative = self.last_matrices
            lines = [
                r"\documentclass[11pt]{article}",
                r"\usepackage[utf8]{inputenc}",
                r"\usepackage[T1]{fontenc}",
                r"\usepackage[spanish]{babel}",
                r"\usepackage{amsmath}",
                r"\usepackage[landscape,margin=1.7cm]{geometry}",
                r"\allowdisplaybreaks",
                r"\begin{document}",
                r"\section*{Cinemática directa: Denavit--Hartenberg modificada}",
                r"\subsection*{Transformaciones consecutivas}",
            ]
            for i, matrix in enumerate(individual, start=1):
                lines.append(rf"\[\scriptstyle {{}}^{{{i-1}}}T_{{{i}}}={sp.latex(matrix)}\]")
            lines.append(r"\subsection*{Transformaciones acumuladas desde 0}")
            for i, matrix in enumerate(cumulative, start=1):
                lines.append(rf"\[\scriptstyle {{}}^{{0}}T_{{{i}}}={sp.latex(matrix)}\]")
            lines.append(r"\end{document}")
            try:
                Path(path).write_text("\n".join(lines), encoding="utf-8")
                self.status.set(f"Matrices LaTeX exportadas a {path}")
            except OSError as ex:
                messagebox.showerror("No se pudo exportar", str(ex))

    root = tk.Tk()
    RobotDHApp(root)
    root.mainloop()


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        self_test()
    else:
        start_gui()
