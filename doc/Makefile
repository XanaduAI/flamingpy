# Minimal makefile for Sphinx documentation
# You can set these variables from the command line.
SPHINXBUILD   = sphinx-build
SPHINXAPIDPC  = sphinx-apidoc
SOURCEDIR     = ..
BUILDDIR      = _build
DEMOBUILDDIR  = tutorials/_out


# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(SPHINXOPTS) .
# the i18n builder cannot share the environment and doctrees with the others
I18NSPHINXOPTS  = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .


# Put it first so that "make" without argument is like "make help".
help:
	@echo "Please use \`make <target>' where <target> is one of:"
	@echo "  clean      to remove current _build/ dir"
	@echo "  html       to make standalone HTML files"
	@echo "  docs       to perform \`clean' and then \`html' together"

.PHONY: help


.PHONY: clean
clean:
	rm -rf $(BUILDDIR)/* $(DEMOBUILDDIR)/*


.PHONY: html
html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished."


.PHONY: docs
docs: clean html
