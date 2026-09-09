TOP_DIR = ../..
include $(TOP_DIR)/tools/Makefile.common

DEPLOY_RUNTIME ?= /kb/runtime
TARGET ?= /kb/deployment

APP_SERVICE = app_service

WRAP_PYTHON_TOOL = wrap_python3
WRAP_PYTHON_SCRIPT = bash $(TOOLS_DIR)/$(WRAP_PYTHON3_TOOL).sh

SRC_PERL = $(wildcard scripts/*.pl)
BIN_PERL = $(addprefix $(BIN_DIR)/,$(basename $(notdir $(SRC_PERL))))
DEPLOY_PERL = $(addprefix $(TARGET)/bin/,$(basename $(notdir $(SRC_PERL))))

SRC_SERVICE_PERL = $(wildcard service-scripts/*.pl)
BIN_SERVICE_PERL = $(addprefix $(BIN_DIR)/,$(basename $(notdir $(SRC_SERVICE_PERL))))
DEPLOY_SERVICE_PERL = $(addprefix $(SERVICE_DIR)/bin/,$(basename $(notdir $(SRC_SERVICE_PERL))))

SRC_SERVICE_PYTHON = $(wildcard service-scripts/*.py)
BIN_SERVICE_PYTHON = $(addprefix $(BIN_DIR)/,$(basename $(notdir $(SRC_SERVICE_PYTHON))))
DEPLOY_SERVICE_PYTHON = $(addprefix $(SERVICE_DIR)/bin/,$(basename $(notdir $(SRC_SERVICE_PYTHON))))

CLIENT_TESTS = $(wildcard t/client-tests/*.t)
SERVER_TESTS = $(wildcard t/server-tests/*.t)
PROD_TESTS = $(wildcard t/prod-tests/*.t)

STARMAN_WORKERS = 8
STARMAN_MAX_REQUESTS = 100

TPAGE_ARGS = --define kb_top=$(TARGET) --define kb_runtime=$(DEPLOY_RUNTIME) --define kb_service_name=$(SERVICE) \
	--define kb_service_port=$(SERVICE_PORT) --define kb_service_dir=$(SERVICE_DIR) \
	--define kb_sphinx_port=$(SPHINX_PORT) --define kb_sphinx_host=$(SPHINX_HOST) \
	--define kb_starman_workers=$(STARMAN_WORKERS) \
	--define kb_starman_max_requests=$(STARMAN_MAX_REQUESTS)

PPI_CONDA_ENV = /opt/conda-cepi-ppi

all: bin 

bin: $(BIN_PERL) $(BIN_SERVICE_PERL) $(BIN_SERVICE_PYTHON)

deploy: deploy-all
deploy-all: deploy-client 
deploy-client: deploy-libs deploy-scripts deploy-docs

deploy-service: deploy-libs deploy-scripts deploy-service-scripts deploy-local-tools deploy-specs

deploy-dir:
	if [ ! -d $(SERVICE_DIR) ] ; then mkdir $(SERVICE_DIR) ; fi
	if [ ! -d $(SERVICE_DIR)/bin ] ; then mkdir $(SERVICE_DIR)/bin ; fi

deploy-docs: 


clean:

$(BIN_DIR)/predict_ppi: service-scripts/predict_ppi.py
	export KB_CONDA_ENV=$(PPI_CONDA_ENV); \
	$(WRAP_PYTHON_SCRIPT) '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@

#
# ppi_report is stdlib-only, but it is wrapped against the same conda env as
# predict_ppi rather than the runtime python: it runs inside the CEPI-PPI
# container as part of the job, where /opt/conda-cepi-ppi is the interpreter
# guaranteed to be present.
#
$(BIN_DIR)/ppi_report: service-scripts/ppi_report.py
	export KB_CONDA_ENV=$(PPI_CONDA_ENV); \
	$(WRAP_PYTHON_SCRIPT) '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@

$(BIN_DIR)/%: service-scripts/%.pl $(TOP_DIR)/user-env.sh
	$(WRAP_PERL_SCRIPT) '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@

$(BIN_DIR)/%: service-scripts/%.py $(TOP_DIR)/user-env.sh
	$(WRAP_PYTHON_SCRIPT) '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@

deploy-local-tools:
	if [ "$(KB_OVERRIDE_TOP)" != "" ] ; then sbase=$(KB_OVERRIDE_TOP) ; else sbase=$(TARGET); fi; \
	export KB_TOP=$(TARGET); \
	export KB_RUNTIME=$(DEPLOY_RUNTIME); \
	export KB_CONDA_ENV=$(PPI_CONDA_ENV); \
	for script in predict_ppi ppi_report ; do \
	    cp service-scripts/$$script.py $(TARGET)/pybin; \
	    $(WRAP_PYTHON_SCRIPT) "$$sbase/pybin/$$script.py" $(TARGET)/bin/$$script; \
	done



include $(TOP_DIR)/tools/Makefile.common.rules
