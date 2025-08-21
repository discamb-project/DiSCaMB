# -*- mode: python ; coding: utf-8 -*-

import vispy.io
import vispy.glsl
import vispy.gloo

data_files = [
    (os.path.dirname(vispy.gloo.__file__), os.path.join("vispy", "gloo")), 
    (os.path.dirname(vispy.glsl.__file__), os.path.join("vispy", "glsl")),
    (os.path.join(os.path.dirname(vispy.io.__file__), "_data"), os.path.join("vispy", "io", "_data"))
]

hidden_imports = [
        #"vispy.ext._bundled.six",
        "vispy.app.backends._pyside6",    
        "vispy.glsl",   
        "vispy.gloo.gl.glplus",
        "vispy.gloo",
    ]

a = Analysis(
    ['dGui2.py'],
    pathex=[],
    binaries=[],
    datas=data_files,
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    exclude_binaries=True,
    name='dGui',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='dGui',
)
