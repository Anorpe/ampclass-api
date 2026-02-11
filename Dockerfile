FROM python:3.9-slim

WORKDIR /multipepgen-app

# Instalar dependencias del sistema para compilar (gcc, etc)
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    && rm -rf /var/lib/apt/lists/*

# Instalar poetry
RUN pip install poetry

# Copiar archivos de dependencias
COPY pyproject.toml ./

# Instalar dependencias (sin crear virtualenv para que queden en el sistema)
RUN poetry config virtualenvs.create false \
    && poetry install --no-root

# Copiar el resto del codigo
COPY . .

# Exponer puerto
EXPOSE 5000

# Comando por defecto
CMD ["python", "api.py"]
