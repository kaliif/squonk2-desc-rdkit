# ---- base python ----
FROM python:3.13.9-slim-trixie as python-base

ENV PYTHONUNBUFFERED=1
ENV HOME=/code


# another stage for poetry installation. this ensures poetry won't end
# up in final image where it's not needed
FROM python-base AS poetry-base

ARG POETRY_VERSION=2.1.4
RUN pip install --no-cache-dir poetry==${POETRY_VERSION}

WORKDIR /
COPY poetry.lock pyproject.toml /

RUN POETRY_VIRTUALENVS_IN_PROJECT=true poetry install --no-root --only main --no-directory



# final stage. only copy the venv with installed packages and point
# paths to it
FROM python-base as final

COPY --from=poetry-base /.venv /.venv

ENV PYTHONPATH="${PYTHONPATH}:/.venv/lib/python3.13/site-packages/"
ENV PATH=/.venv/bin:$PATH


WORKDIR ${HOME}
COPY src/ ./
