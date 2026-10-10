"""Isotherm registry: find isotherms by the name used in project files.

Built-in isotherms register themselves when ``chromsim.isotherms`` is imported. A new isotherm is
a new file with a class decorated with ``@register_isotherm`` (CLAUDE.md rule 6). The registry is
filled at import time only and not changed afterwards.
"""

from chromsim.errors import IsothermConfigError, PluginRegistrationError
from chromsim.isotherms.base import IsothermBase

_ISOTHERMS: dict[str, IsothermBase] = {}


def register_isotherm(isotherm_class: type[IsothermBase]) -> type[IsothermBase]:
    """Class decorator that adds an isotherm to the registry under ``isotherm_class.name``.

    Raises
    ------
    PluginRegistrationError
        If another isotherm already uses the same name.
    """
    name = isotherm_class.name
    if name in _ISOTHERMS:
        raise PluginRegistrationError(
            f"An isotherm named '{name}' is already registered "
            f"({type(_ISOTHERMS[name]).__qualname__}); choose a different name."
        )
    _ISOTHERMS[name] = isotherm_class()
    return isotherm_class


def get_isotherm(name: str) -> IsothermBase:
    """Return the isotherm registered under ``name``.

    Raises
    ------
    IsothermConfigError
        If no isotherm has that name; the message lists the available ones.
    """
    if name not in _ISOTHERMS:
        raise IsothermConfigError(
            f"unknown isotherm '{name}'. Available isotherms: {', '.join(list_isotherms())}."
        )
    return _ISOTHERMS[name]


def list_isotherms() -> list[str]:
    """Names of all registered isotherms, sorted."""
    return sorted(_ISOTHERMS)
