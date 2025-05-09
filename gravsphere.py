# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_gravsphere')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_gravsphere')
    _gravsphere = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_gravsphere', [dirname(__file__)])
        except ImportError:
            import _gravsphere
            return _gravsphere
        try:
            _mod = imp.load_module('_gravsphere', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _gravsphere = swig_import_helper()
    del swig_import_helper
else:
    import _gravsphere
del _swig_python_version_info

try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except __builtin__.Exception:
    class _object:
        pass
    _newclass = 0


def GetLikeBinAnisSphVSPfunc(R_list, sigmas, sigerrs, rho_params, beta_params, plum_params, vir_shape, stellar_mass_3r, rh, rows, min_rad, max_rad):
    return _gravsphere.GetLikeBinAnisSphVSPfunc(R_list, sigmas, sigerrs, rho_params, beta_params, plum_params, vir_shape, stellar_mass_3r, rh, rows, min_rad, max_rad)
GetLikeBinAnisSphVSPfunc = _gravsphere.GetLikeBinAnisSphVSPfunc

def GetLikeBinAnisSphVSPPowerLawfunc(R_list, sigmas, sigerrs, rho_params, beta_params, plum_params, vir_shape, stellar_mass_3r, rh, rows, min_rad, max_rad):
    return _gravsphere.GetLikeBinAnisSphVSPPowerLawfunc(R_list, sigmas, sigerrs, rho_params, beta_params, plum_params, vir_shape, stellar_mass_3r, rh, rows, min_rad, max_rad)
GetLikeBinAnisSphVSPPowerLawfunc = _gravsphere.GetLikeBinAnisSphVSPPowerLawfunc

def PowerLawFitfunc(R_list, sigmas, sigerrs, rho_params, beta_params, plum_params, vir_shape, stellar_mass_3r, rh, rows, results, vsp1, vsp2, min_rad, max_rad):
    return _gravsphere.PowerLawFitfunc(R_list, sigmas, sigerrs, rho_params, beta_params, plum_params, vir_shape, stellar_mass_3r, rh, rows, results, vsp1, vsp2, min_rad, max_rad)
PowerLawFitfunc = _gravsphere.PowerLawFitfunc

def ZhaoFitfunc(R_list, sigmas, sigerrs, rho_params, beta_params, plum_params, vir_shape, stellar_mass_3r, rh, rows, results, vsp1, vsp2, min_rad, max_rad):
    return _gravsphere.ZhaoFitfunc(R_list, sigmas, sigerrs, rho_params, beta_params, plum_params, vir_shape, stellar_mass_3r, rh, rows, results, vsp1, vsp2, min_rad, max_rad)
ZhaoFitfunc = _gravsphere.ZhaoFitfunc
# This file is compatible with both classic and new-style classes.


