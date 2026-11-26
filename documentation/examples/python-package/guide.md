---
title: Read a saved Workbench request
author: Gunwoo Park
version: '0.2.0'
date: '2026-11-26'
lang: en
---

[AI-generated text](#ai-assistance)

## Purpose {#purpose}

Your research script has the receipt from an earlier Workbench request. Use
`wb-python` to attach to the selected existing Thread from an ordinary Python
process, read that request's status, and retrieve its result when complete.
Keeping the original receipt lets the script identify the exact request on
subsequent runs.

A Session Profile groups a workspace's work and recovery records. A Thread
identifies one conversation within that Profile. A receipt identifies the
specific request inside the Thread.

## Installation {#installation}

Read or render this guide without a Workbench runtime. To run the example,
use Python 3.11 or later in a configured Sandbox and a writable directory.
Set `WB_PYTHON_SOURCE` to a complete wb-python source directory or wheel
supplied by the operator. In a full development checkout you can use its
`wb-python/` source package. The public `wb-python/` entry contains the manual
and references. With pip available, install the runtime into a fresh
virtual environment:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install "$WB_PYTHON_SOURCE"
```

Set the following environment variables from the existing operator and
request records. Keep the actual values in your own environment.

| Setting | Value to supply |
| --- | --- |
| `WB_WORKSPACE`, `WB_HELPER`, `WB_HELPER_SHA256` | Configured workspace root, installed Session Protocols helper path, and its reviewed SHA-256 pin |
| `WB_PROFILE`, `WB_THREAD` | Selected Session Profile name and the exact admitted Thread already confirmed for it |
| `WB_SOCKET` | Existing local socket for an operator-selected supported native runtime (Codex 0.171.0 or 0.172.0) |
| `WB_RECEIPTS`, `WB_RECEIPT` | The original private SDK receipt directory and the 64-character receipt hash saved with the request |

## First useful run {#first-use}

Discover the existing Profile names before selecting `WB_PROFILE` and its
confirmed `WB_THREAD`:

```sh
wb-python --workspace "$WB_WORKSPACE" --helper "$WB_HELPER" \
  --helper-sha256 "$WB_HELPER_SHA256" profiles
```

Save the following as `read_retained_request.py` and run
`python read_retained_request.py` in the same environment.

```python
import json
import os

from wb_python import ProfileAccess, ReceiptLog, Workbench

profile = os.environ["WB_PROFILE"]
thread = os.environ["WB_THREAD"]
profiles = ProfileAccess(
    os.environ["WB_WORKSPACE"], os.environ["WB_HELPER"],
    os.environ["WB_HELPER_SHA256"],
)
reference = ReceiptLog(os.environ["WB_RECEIPTS"]).read(os.environ["WB_RECEIPT"])

with Workbench(profiles, os.environ["WB_SOCKET"], allowed_threads=(thread,),
               enabled_operations=("read",)) as workbench:
    selected = workbench.attach(profile, thread, resume=False)
    if reference.session != selected:
        raise SystemExit("Receipt belongs to a different Session reference.")

    snapshot = workbench.status(reference)
    summary = {"state": snapshot["state"]}
    if snapshot["state"] == "completed":
        result = workbench.result(reference)
        summary.update(message_count=len(result["messages"]),
                       code_execution=result["code_execution"])
    print(json.dumps(summary))
```

For a completed request with one response, the script prints
`{"state": "completed", "message_count": 1, "code_execution": "inert_text_only"}`.
The message count depends on the original request. `result["messages"]` holds
the text for the calling script to review. The example prints only the count
and uses native read operations.

## API reference {#api-reference}

| Callable | Input | Result |
| --- | --- | --- |
| `ProfileAccess(workspace, helper_path, helper_sha256)` | Operator-selected configuration and verified helper pin | Access to existing Session Profiles |
| `ReceiptLog(directory).read(receipt_hash)` | Original private SDK receipt store and hash | The retained request reference |
| `Workbench.attach(profile, thread_id, resume=False)` | Chosen existing Profile and admitted Thread | Validated Session reference |
| `Workbench.status(reference)` and `Workbench.result(reference)` | Exact retained request reference | Current request state and its terminal result when available |

## Troubleshooting and scope {#troubleshooting}

If the state is `queued`, `inProgress`, or `unresolved`, retain the original
receipt and check it later. `failed` and `interrupted` identify terminal
outcomes. If the binding or receipt differs, confirm the intended Profile,
Thread and receipt with the operator before running again. The rendered
example documents this configured workflow. It contains no workspace or
receipt data.

## Author and project information {#author-and-project-information}

Gunwoo Park directs Agentic Workbench and its `wb-templates` collection.

Project overview: <https://github.com/gwpark-git/agentic-workbench-knowledge>  
wb-templates repository (planned): <https://github.com/gwpark-git/wb-templates>

### AI assistance

OpenAI Codex assisted with drafting this example under Gunwoo Park's direction.
