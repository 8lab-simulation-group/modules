## $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/modules/module-aerodyn/Makefile.inc,v 1.5 2009/06/12 11:14:56 masarati Exp $
MODULE_DEPENDENCIES= mbdyn2AD.lo 
MODULE_LINK=-L. -lAeroDyn
# FIXME: needed to link objects compiled with Intel ifc
#MODULE_LINK=-L/opt/intel/ia32/lib/ -lF90 -lIEPCF90 -lintrins $(MODULE_LINK)
