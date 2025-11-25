
# 🚀 BLinCS Developer Setup Guide  
*A complete guide for developers working with BLinCS — reliable, step‑by‑step, cross‑platform.*

This document explains **exactly** how to set up BLinCS as a developer.  
It handles common issues such as wrong Python interpreter, broken virtual envs, missing MkDocs, etc.

---

## 1. 📦 Clone the repository

Open a terminal (macOS/Linux) or PowerShell (Windows):

```bash
git clone https://github.com/turbulencelover/BLinCS.git
cd BLinCS
```

You should now see:

```
BLinCS/
    src/
    docs/
    pyproject.toml
    ...
```

---

## 2. 🐍 Create a clean Python virtual environment

Always create a project‑local venv. Never rely on system Python or Anaconda.

### macOS / Linux
```bash
python3 -m venv .venv
```

### Windows (PowerShell)
```powershell
python -m venv .venv
```

Confirm it exists:

```bash
ls .venv
```

---

## 3. 🧠 Two safe ways to use the virtual environment

You can **activate** it or **use full paths**.

### ✔ Method A — Activate the venv

macOS/Linux:

```bash
source .venv/bin/activate
```

Windows:

```powershell
.\.venv\Scripts\activate
```

Prompt should show `(.venv)`.

### ✔ Method B — Recommended (avoids activation issues)

Use full paths directly:

```bash
./.venv/bin/python
./.venv/bin/mkdocs
```

Example:

```bash
./.venv/bin/python -c "print('venv working')"
```

This is guaranteed to use the correct environment.

---

## 4. 📥 Install BLinCS in editable developer mode

```bash
./.venv/bin/python -m pip install -e .
```

Test the installation:

```bash
./.venv/bin/python -c "import blincs; print(blincs.__file__)"
```

Expected output:

```
BLinCS/src/blincs/__init__.py
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

If this runs, your installation is correct.

---

## 6. 🧪 Running tests

Tests are located under:

```
tests/
```

Run a test:

```bash
./.venv/bin/python tests/NonlinearEffects/<test_script>.py
```

---

## 7. 📝 Editing source code (developer workflow)

All source code is under:

```
src/blincs/
```

Because BLinCS was installed with `pip install -e .`, code edits take effect immediately — **no reinstall needed**.

---

## 8. 📚 Build and preview documentation (MkDocs)

### Install doc dependencies:

```bash
./.venv/bin/python -m pip install -r docs/requirements.txt
```

### Serve the docs locally:

```bash
./.venv/bin/mkdocs serve
```

Open:

👉 http://127.0.0.1:8000

Live reload works automatically.

---

## 9. 🧹 Common troubleshooting

### ❌ `zsh: command not found: python`
Use:

```bash
python3
```
or:

```bash
./.venv/bin/python3
```

---

### ❌ `mkdocs: command not found`
You installed MkDocs globally. Install it inside the venv:

```bash
./.venv/bin/python -m pip install -r docs/requirements.txt
```

---

### ❌ `ModuleNotFoundError: No module named "blincs"`
You forgot to install the package:

```bash
./.venv/bin/python -m pip install -e .
```

---

### ❌ Wrong interpreter used in VS Code or Spyder

macOS/Linux:

```
/path/to/BLinCS/.venv/bin/python3
```

Windows:

```
BLinCS\.venv\Scripts\python.exe
```

---

### ❌ Broken venv (common after uninstalling Anaconda)

Just recreate it:

```bash
rm -rf .venv
python3 -m venv .venv
./.venv/bin/python -m pip install -e .
```

---

## 10. 🌐 IDE Setup (VS Code & Spyder)

### VS Code

1. Ctrl+Shift+P → **Python: Select Interpreter**
2. Choose:
   ```
   BLinCS/.venv/bin/python3
   ```

### Spyder

1. Tools → Preferences → Python Interpreter
2. Select:
   ```
   BLinCS/.venv/bin/python
   ```
3. Install Spyder kernels in your venv:
   ```bash
   ./.venv/bin/python -m pip install "spyder-kernels>=3.0,<3.1"
   ```
4. Restart kernel.

---

## 11. 🔄 Keeping your copy up to date

```bash
git pull
```

If dependencies changed:

```bash
./.venv/bin/python -m pip install -e .
./.venv/bin/python -m pip install -r docs/requirements.txt
```

---

## 12. 🎉 You're ready to develop BLinCS!

You can now:

- run simulations  
- edit/extend the code  
- build documentation  
- run tests  
- contribute pull requests  

Welcome to BLinCS development 🚀
