# wb-python interface preset

[AI-generated text](#ai-assistance)

[Complete wb-python manual and reference](../wb-python/README.md)

The `wb-python` preset gives Python scripts access to selected Workbench
operations. It provides a reusable integration base for future
`wb-applications` packages that need to work with profiles and native threads
without building a new session registry.

| Script task | Preset interface | Scope |
|---|---|---|
| Inspect Session Profiles | `ProfileAccess.list()` uses the existing workspace's approved Profile helper. | Read metadata for authorized profiles. Access to the helper and workspace comes from the local operator. |
| Reuse an existing conversation | `Workbench.attach(..., resume=False)` checks a selected profile and native thread reference. | Attachment can be read-only. Explicit resume and request submission use the configured local runtime. |
| Select an Agent Profile for a task | `Workbench.create_session(profile, agent_profile, task)` uses the configured local bootstrap. | Use an explicitly selected existing Agent Profile. This operation can start a native session and requires the configured bootstrap and the operator's authorization. |

A **Session Profile** groups work and recovery records in a workspace.
An **Agent Profile** supplies the persona and invocation settings for an
agent. Both stay under the existing Workbench lifecycle. Returned code and
messages are available for inspection as text.

The maintained preset lives in `wb-templates` and is planned for an independent
release. This collection includes the manual, API reference and a worked
readback guide. For execution, obtain the full runtime separately and follow
the Sandbox setup in the [Reader Guide](../wb-python/docs/READER_GUIDE.org). `wb-client`
is a separate project for Host-side access and its own release.

## Jupyter integration in wb-applications

`wb-applications` owns the separate `wb-jupyter` package. It carries its own
reader's guide, dependency declarations, tests and clean first-use notebook.
Select it when a workflow needs the IPython `%wb` and `%%wb` commands.
These commands use a configured `wb-client` to select an existing Profile
and Thread, submit explicit requests, and inspect the result. Inside the
Sandbox, the configured Client can use `wb-python`'s local transport.

The notebook starts in an unconfigured mode. Reviewing returned text and
inserting an unexecuted cell are distinct actions in the supported
Lab/Kernel configuration.

A standalone release of `wb-jupyter` is under consideration. For most
Host-facing wb-applications packages, the release plan is individual
distribution after security review.

[Read a retained Workbench request from a Python script](examples/python-package/preview.html)

[Python package guide source and other documentation formats](README.md#documentation-families)

[Collection guide and functional reference](collection-guide/report.html)

## Author and project information

**Author and Project Director:** Gunwoo Park.

Project overview: <https://github.com/gwpark-git/agentic-workbench-knowledge>  
wb-templates repository (planned): <https://github.com/gwpark-git/wb-templates>

### AI assistance

OpenAI Codex assisted with drafting and editing this documentation under
Gunwoo Park's direction.
