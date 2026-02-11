FROM python:3.12-slim

# Maintainer info
LABEL maintainer="aorregope@gmail.com"

# Make working directories
RUN mkdir -p /multipepgen-app
WORKDIR /multipepgen-app

# Install poetry
RUN pip install poetry

# Copy poetry files
COPY pyproject.toml poetry.lock ./

# Install dependencies
RUN poetry install --no-root

# Copy the rest of the application code
COPY . .

# Run the python application
#CMD ["poetry", "run", "gunicorn", "-w", "4", "-b", "0.0.0.0:8050", "app:server"]
CMD ["poetry", "run", "gunicorn","-b", "0.0.0.0:8050", "app:server"]
