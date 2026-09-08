"""build_core() for viva-tumor-tcell.

Builds on viva-munk's fully-registered core (pymunk_agent types + PymunkProcess
+ the `_add`/`_remove` realize and tuple-delta apply patches) and adds this
workspace's own ``tumor_tcell_agent`` type and processes.
"""
from viva_munk.core import core_import

from .types import register_types
from .processes.tumor import TumorCellProcess
from .processes.t_cell import TCellProcess
from .processes.dendritic_cell import DendriticCellProcess
from .processes.physics import TumorTcellPhysics
from .processes.field import DiffusionField


def build_core(core=None):
    # viva-munk's core registers pymunk_agent, PymunkProcess, the set_float /
    # concentration types, and the _add/realize + tuple-delta apply patches.
    core = core_import(core)
    register_types(core)
    core.register_link('TumorCellProcess', TumorCellProcess)
    core.register_link('TCellProcess', TCellProcess)
    core.register_link('DendriticCellProcess', DendriticCellProcess)
    core.register_link('TumorTcellPhysics', TumorTcellPhysics)
    core.register_link('DiffusionField', DiffusionField)
    return core
