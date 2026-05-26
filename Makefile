# =============================================================================
# Makefile - Automatisation du projet hybride Rust/Python (Maturin)
# =============================================================================

# Détection de Python (cherchera 3.14 en priorité, puis 3, puis python)
PYTHON_CMD ?= $(shell command -v python3.14 || command -v python3 || command -v python)

API_DIR     := python
SANDBOX_DIR := sandbox_python
VENV_DIR    := $(SANDBOX_DIR)/.venv

# Chemins absolus vitaux pour ne pas se faire "hijacker" par le shell parent
VENV_ABS     := $(abspath $(VENV_DIR))
VENV_PYTHON  := $(VENV_ABS)/bin/python
VENV_PIP     := $(VENV_ABS)/bin/pip
VENV_MATURIN := $(VENV_ABS)/bin/maturin

.PHONY: all help build clean rebuild venv run

all: build

# 1. Création de l'environnement virtuel
venv:
	@echo "[VENV] Création de l'environnement virtuel dans $(VENV_DIR)..."
	@$(PYTHON_CMD) -m venv $(VENV_DIR)
	@echo "[VENV] Mise à jour de pip et installation des outils..."
	@$(VENV_PIP) install --upgrade pip maturin
	@echo "[VENV] Installation des dépendances de $(SANDBOX_DIR)..."
	@cd $(SANDBOX_DIR) && $(VENV_PIP) install -e . || echo "[WARN] Aucun fichier d'installation trouvé."

# 2. Compilation blindée (Contournement de la détection capricieuse de Maturin)
build: venv
	@echo "[BUILD] Nettoyage des anciens artefacts de build..."
	@rm -rf $(API_DIR)/target/wheels
	@echo "[BUILD] Compilation de l'API avec ciblage strict..."
	@# L'argument '-i' force Maturin à compiler UNIQUEMENT pour notre interpréteur cible
	@cd $(API_DIR) && $(VENV_MATURIN) build --release -i $(VENV_PYTHON) --out target/wheels
	@echo "[BUILD] Injection de la librairie compilée dans le VENV local..."
	@# On force pip à installer le wheel fraîchement généré pour écraser tout conflit
	@$(VENV_PIP) install --force-reinstall --no-deps $(API_DIR)/target/wheels/ferrispline*.whl
	@echo "[BUILD] Compilation et installation terminées avec succès."

# 3. Exécution
re: build
	@echo "[RUN] Exécution du programme principal..."
	@$(VENV_PYTHON) $(SANDBOX_DIR)/src/nurbs_math/main.py


# Permet de passer le fichier et optionnellement le nombre de samples : make run mon_fichier.vtk 200
ifeq (run,$(firstword $(MAKECMDGOALS)))
  FILE_ARG := $(word 2,$(MAKECMDGOALS))
  SAMPLES_ARG := $(word 3,$(MAKECMDGOALS))
  $(eval $(FILE_ARG):;@:)
  $(eval $(SAMPLES_ARG):;@:)
endif

run:
	@if [ -z "$(FILE_ARG)" ]; then \
		echo "ERROR: FILE is not set. Usage: make run path/to/file.vtk [samples]"; \
		exit 1; \
	fi
	@echo "[RUN] Exécution avec le fichier : $(FILE_ARG) et samples : $(if $(SAMPLES_ARG),$(SAMPLES_ARG),100)"
	@$(VENV_PYTHON) $(SANDBOX_DIR)/src/nurbs_math/main.py --file $(FILE_ARG) --samples $(if $(SAMPLES_ARG),$(SAMPLES_ARG),100)

# 4. Nettoyage
clean:
	@echo "[CLEAN] Suppression de l'environnement virtuel et des artefacts..."
	@rm -rf $(VENV_DIR)
	@rm -rf $(API_DIR)/target
	@cargo clean --manifest-path $(API_DIR)/Cargo.toml 2>/dev/null || true
	@cargo clean --manifest-path core_rust/Cargo.toml 2>/dev/null || true
	@find . -type d -name "__pycache__" -exec rm -rf {} +
	@find . -type f -name "*.pyc" -delete
	@find . -type f -name "*.so" -delete
	@find . -type f -name "*.pyd" -delete
	@find . -type f -name "*.whl" -delete
	@echo "[CLEAN] Environnement nettoyé."

rebuild: clean build

help:
	@echo "Commandes disponibles :"
	@echo "  make venv    : Crée l'environnement virtuel."
	@echo "  make build   : Compile l'API via Maturin (ciblage strict) et l'installe."
	@echo "  make run     : Compile (si nécessaire) puis lance l'application."
	@echo "  make clean   : Détruit le venv et les caches."
	@echo "  make rebuild : Nettoyage complet puis compilation."
