"""Optimizer registry: find optimizers by name. A new optimizer is a new file with a class
decorated with ``@register_optimizer`` (see docs/how_to/add_optimizer.md). Filled at import time
only."""

from chromsim.errors import ChromSimError, PluginRegistrationError

_OPTIMIZERS: dict[str, type] = {}


def register_optimizer(optimizer_class: type) -> type:
    """Class decorator adding an optimizer under ``optimizer_class.name``."""
    name = getattr(optimizer_class, "name")  # noqa: B009  (any class with a `name` attribute)
    if name in _OPTIMIZERS:
        raise PluginRegistrationError(f"An optimizer named '{name}' is already registered.")
    _OPTIMIZERS[name] = optimizer_class
    return optimizer_class


def get_optimizer(name: str, **options):
    """A new instance of the optimizer ``name`` with the given options."""
    if name not in _OPTIMIZERS:
        raise ChromSimError(
            f"Unknown optimizer '{name}'. Available: {', '.join(list_optimizers())}."
        )
    return _OPTIMIZERS[name](**options)


def list_optimizers() -> list[str]:
    """Names of all registered optimizers, sorted."""
    return sorted(_OPTIMIZERS)
