# MultiPepGen

## Docker

### Build

```docker build -t multipepgen-app-image .```
### Run

```docker run --network="host" -p 8050:8050 multipepgen-app-image```   
### Test

```bash curl -X POST "http://localhost:8050/multipepgen" -H "accept: application/json" -H "Content-Type: application/json" -d "{ \"peptides\": [ \"PEPTIDE1\", \"PEPTIDE2\" ]}"```

## Command line

```poetry run gunicorn -b 0.0.0.0:8050 app:server   ```

## Local
### Install

```poetry install``` 
### Run

```poetry run python app.py```

### TMUX
Para el manejo de sesiones de terminal en segundo plano se utiliza tmux.

crear sesión
```tmux new -s multipepgen```

Salir de la sesion
```Ctrl + b + d```
o
```tmux detach```

ver sesiones activas
```tmux ls```

volver a retomar una sesion
```tmux attach -t multipepgen```
o
```tmux a -t multipepgen```

Eliminar una sesion
```tmux kill-session -t multipepgen``` 
