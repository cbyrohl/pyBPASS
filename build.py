from numpy.distutils.core import Extension, setup

def build(setup_kwargs):
    sed_ext = Extension(
                    name="pyBPASS._sed_tools",
                    sources=["fsrc/sed_tools.f90"]
    )
    ext_modules = [sed_ext]
    setup_kwargs.update({
        "ext_modules": ext_modules,
        "zip_safe": False,
    })
