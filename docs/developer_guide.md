
# 🚀 BLinCS Developer Setup Guide  
*A complete guide for developers working with BLinCS — reliable, step‑by‑step, cross‑platform.*

This document explains **exactly** how to set up BLinCS as a developer.
It handles common issues such as wrong Python interpreter, broken virtual envs, missing MkDocs, etc.

---

## 1. 📦 Clone the repository

```bash
git clone https://github.com/turbulencelover/BLinCS.git
cd BLinCS
```

---

## 2. 🐍 Create a clean Python virtual environment

### macOS / Linux
```bash
python3 -m venv .venv
```

### Windows (PowerShell)
```powershell
python -m venv .venv
```

Check:

```bash
ls .venv
```

---

## 3. 🧠 Two ways to use the venv

### Method A — Activate it

macOS/Linux:
```bash
source .venv/bin/activate
```

Windows:
```powershell
.\.venv\Scripts\activate
```

### Method B — Recommended: Use full paths

```bash
./.venv/bin/python
./.venv/bin/mkdocs
```

This avoids activation issues.

---

## ⚠️ 3.1 Important Note About Virtual Environments

A new or recreated `.venv` starts empty — no BLinCS, no MkDocs.
If VS Code or Spyder switches interpreters, you may suddenly lose access to `blincs`.

### ✔ After creating or changing `.venv`, you MUST reinstall BLinCS:

```bash
./.venv/bin/python -m pip install -e .
```

Or in VS Code:

1. Select interpreter: `BLinCS/.venv/bin/python`
2. Command Palette → Python: Create Terminal
3. Then run:
   ```bash
   python -m pip install -e .
   ```

If you skip this, you'll get:
```
ModuleNotFoundError: No module named 'blincs'
```

### ✔ The interpreter that runs your code must be the interpreter where you installed BLinCS.

Check:

```bash
python -c "import sys; print(sys.executable)"
```

Must show something like:

```
BLinCS/.venv/bin/python
```

---

## 4. 📥 Install BLinCS in editable mode

```bash
./.venv/bin/python -m pip install -e .
```

Test:

```bash
./.venv/bin/python -c "import blincs; print(blincs.__file__)"
```

---

## 5. ▶️ Run a quick simulation test

```bash
./.venv/bin/python - << 'EOF'
from blincs.blwave_simulation import simulate_BLwave
from blincs.vel_base import LogProfile
from blincs.forcing import BoxFarmGeom, TurbineBoxForcing

geom = BoxFarmGeom(ax=20000, az=100, zh=100)
forcing = TurbineBoxForcing(geom)
bg = LogProfile(U_top=20, z0=1, ztop=1000)

x,z,ur,wr,eta,fx,u0z = simulate_BLwave(
    g_prime=0.3, N=0.01,
    background=bg, forcing=forcing,
    Nx=200, Nz=100, Lx=1e6, Lz=1000
)
print("BLinCS simulation ran successfully!")
EOF
```

---

## 6. 🧪 Running tests

```bash
./.venv/bin/python tests/NonlinearEffects/<test>.py
```

---

## 7. 📝 Editing BLinCS source

Your code lives in:

```
src/blincs/
```

Because you installed in editable mode, changes take effect immediately.

---

## 8. 📚 Build documentation (MkDocs)

Install doc deps:

```bash
./.venv/bin/python -m pip install -r docs/requirements.txt
```

Serve docs:

```bash
./.venv/bin/mkdocs serve
```

Visit:

```
http://127.0.0.1:8000
```

---

## 9. 🧹 Troubleshooting

### ❌ `python: command not found`
Use:
```bash
python3
```
or:
```bash
./.venv/bin/python
```

### ❌ `mkdocs: command not found`
Install in venv:
```bash
./.venv/bin/python -m pip install -r docs/requirements.txt
```

### ❌ `ModuleNotFoundError: blincs`
Reinstall:
```bash
./.venv/bin/python -m pip install -e .
```

### ❌ Using the wrong interpreter
Must be:
```
BLinCS/.venv/bin/python
```

### ❌ Broken venv
Recreate:

```bash
rm -rf .venv
python3 -m venv .venv
./.venv/bin/python -m pip install -e .
```

---

## 10. 🌐 IDE Setup

### VS Code
1. Ctrl+Shift+P → Python: Select Interpreter
2. Choose:
   ```
   BLinCS/.venv/bin/python
   ```

### Spyder
1. Preferences → Python Interpreter
2. Select:
   ```
   BLinCS/.venv/bin/python
   ```
3. Install kernel:
   ```bash
   ./.venv/bin/python -m pip install "spyder-kernels>=3.0,<3.1"
   ```

---

## 11. 🔄 Updating BLinCS

```bash
git pull
```

If deps changed:

```bash
./.venv/bin/python -m pip install -e .
./.venv/bin/python -m pip install -r docs/requirements.txt
```

---

## 12. 🎉 Ready to Develop BLinCS!

You can now:

- run simulations  
- update source code  
- build documentation  
- run tests  
- contribute pull requests  

Welcome to BLinCS development 🚀
