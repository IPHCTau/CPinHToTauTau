import os
import law
import order as od
from typing import Optional

def add_expanded_shift_aliases(
        config: od.Config,
        shift_source: str,
        aliases: dict,
        variations: Optional[list]=['up','down'],
) -> None:
    for direction in variations:
        shift = config.get_shift(od.Shift.join_name(shift_source, direction))
        _aliases = shift.x("column_aliases", {})
        # format keys and values
        inject_shift = lambda s: re.sub(r"\{([^_])", r"{_\1", s).format(**shift.__dict__)
        _aliases.update({inject_shift(key): inject_shift(value) for key, value in aliases.items()})
        # extend existing or register new column aliases
        shift.x.column_aliases = _aliases
