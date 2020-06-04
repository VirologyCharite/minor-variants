.PHONY: check, tcheck, flake8, clean

XARGS := xargs $(shell test $$(uname) = Linux && echo -r)
DATE := $(shell date +%Y-%m-%d)

check:
	env PYTHONPATH=. python -m discover -v

tcheck:
	env PYTHONPATH=. trial --rterrors test

flake8:
	find bin code test -name '*.py' -print0 | $(XARGS) -0 flake8 --ignore E402,W504

clean:
	rm -fr _trial_temp
