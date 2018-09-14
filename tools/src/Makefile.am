SUBDIRS = modules custom

bin_SCRIPTS = polaris-run polaris-plot polaris-gen polaris-remote polaris-extra polaris-test
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = polaris-run.in polaris-plot.in polaris-gen.in polaris-remote.in polaris-extra.in polaris-test.in

do_substitution = sed -e 's,[@]POLARIS[@],$(prefix)/,g'

polaris-run: polaris-run.in Makefile
	$(do_substitution) < $(srcdir)/polaris-run.in > polaris-run
	chmod +x polaris-run

polaris-plot: polaris-plot.in Makefile
	$(do_substitution) < $(srcdir)/polaris-plot.in > polaris-plot
	chmod +x polaris-plot

polaris-gen: polaris-gen.in Makefile
	$(do_substitution) < $(srcdir)/polaris-gen.in > polaris-gen
	chmod +x polaris-gen

polaris-remote: polaris-remote.in Makefile
	$(do_substitution) < $(srcdir)/polaris-remote.in > polaris-remote
	chmod +x polaris-remote

polaris-extra: polaris-extra.in Makefile
	$(do_substitution) < $(srcdir)/polaris-extra.in > polaris-extra
	chmod +x polaris-extra

polaris-test: polaris-test.in Makefile
	$(do_substitution) < $(srcdir)/polaris-test.in > polaris-test
	chmod +x polaris-test
