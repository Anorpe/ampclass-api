# MultiPepGen

## Docker

### Build

```docker build -t multipepgen-app-image .```
### Run

```docker run --name multipepgen-app-container --network="host" -p 8100:8100 multipepgen-app-image```   


Ejecutar sin bloquear la terminal

```docker run --name multipepgen-app-container --network="host"  -p 8100:8100 -d multipepgen-app-image```

Detener la ejecucion del contenedor

```docker ps```

```docker stop multipepgen-app-container```


### Test

```bash curl -X POST "http://localhost:8050/multipepgen" -H "accept: application/json" -H "Content-Type: application/json" -d "{ \"peptides\": [ \"PEPTIDE1\", \"PEPTIDE2\" ]}"```

## Command line

```poetry run gunicorn -b 0.0.0.0:8100 app:server   ```

## Local
### Install

```poetry install``` 
### Run

```poetry run python app.py```

### TMUX
Para el manejo de sesiones de terminal en segundo plano se utiliza tmux.

crear sesión
```tmux new -s multipepgen-app-session```

Salir de la sesion
```Ctrl + b + d```
o
```tmux detach```

ver sesiones activas
```tmux ls```

volver a retomar una sesion
```tmux attach -t multipepgen-app-session```
o
```tmux a -t multipepgen-app-session```

Eliminar una sesion
```tmux kill-session -t multipepgen-app-session``` 
