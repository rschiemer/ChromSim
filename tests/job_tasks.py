"""Tiny job tasks for runner tests (module-level, so worker processes can import them)."""

import time


def quick_task(context, name):
    context.report(0.5, "half way")
    context.storage.save_record("outputs", name, {"ok": True})
    return name


def failing_task(context):
    raise ValueError("deliberate failure")


def cancellable_task(context, seconds):
    deadline = time.monotonic() + seconds
    step = 0
    while time.monotonic() < deadline:
        step += 1
        context.report(min(step / 100, 0.99), f"step {step}")  # raises once a cancel is requested
        time.sleep(0.05)
    return "not cancelled"
