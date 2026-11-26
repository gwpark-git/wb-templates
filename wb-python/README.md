# wb-python

[AI-generated text](#ai-assistance)

Use wb-python to bring selected Workbench operations into an ordinary Python
script. A research script can select an existing conversation, submit an
explicit request when authorized, save its acknowledgement, and return later
to read the same request's result.

- A **Session Profile** groups persistent work and recovery records. A
  **Thread** identifies the conversation the script selects.
- An **Agent Profile** describes an agent's reusable role and invocation
  procedure. The local Python API can ask an already configured bootstrap
  to create a child Thread under an existing Session Profile.
- A **receipt** identifies an acknowledged request after the submitting
  process exits. Use it to read status and the selected result.

Start with the [Reader Guide](docs/READER_GUIDE.org). It walks through setup,
read-only receipt recovery, explicit submission, and recovery after uncertain
acceptance. The [API reference](docs/API.org) covers native identities and
method behavior. See also the [operating and data boundaries](docs/SECURITY.org),
[restricted endpoint](docs/REMOTE.org), and [validation scope](docs/VALIDATION.org).

wb-python runs inside a configured Sandbox on POSIX with Python 3.11 or later.
The independent wb-client package is the optional Host-side client. The
separate wb-jupyter package adds IPython commands and uses the Client.
The documentation in this collection covers wb-python 0.2.0. Installation
uses a separately obtained, reviewed copy of the wb-python runtime. A standalone
release is planned. The independently supplied runtime carries its own MIT grant.

See the [collection documentation](../documentation/README.md),
[MIT license](../LICENSE) and [license scope](../LICENSE_SCOPE.md).

## Author and project information

Gunwoo Park directs the [Agentic Workbench project](https://github.com/gwpark-git/agentic-workbench-knowledge)
and its wb-python preset.

### AI assistance

OpenAI Codex helped select, edit and organize the public reader guide and
references under Gunwoo Park's direction.
