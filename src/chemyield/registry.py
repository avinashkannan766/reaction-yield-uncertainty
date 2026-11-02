"""Tiny name -> class registry so configs can pick components by name."""

from __future__ import annotations


class Registry:
    def __init__(self, kind: str):
        self.kind = kind
        self._items: dict[str, type] = {}

    def register(self, name: str):
        def deco(cls):
            self._items[name] = cls
            cls.registry_name = name
            return cls

        return deco

    def get(self, name: str):
        if name not in self._items:
            raise KeyError(f"Unknown {self.kind} '{name}'. Available: {sorted(self._items)}")
        return self._items[name]

    def build(self, spec: dict, **extra):
        spec = dict(spec)
        name = spec.pop("name")
        return self.get(name)(**spec, **extra)

    def names(self) -> list[str]:
        return sorted(self._items)
