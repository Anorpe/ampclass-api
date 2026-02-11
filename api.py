from flask import Flask, request, jsonify
import predictor
import pandas as pd
import os

app = Flask(__name__)

# Configurar CORS manualmente para evitar problemas de dependencias
@app.after_request
def after_request(response):
    response.headers.add('Access-Control-Allow-Origin', '*')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
    response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE,OPTIONS')
    return response

@app.route('/predict', methods=['POST'])
def predict():
    try:
        data = request.get_json()
        sequences = data.get('sequences', [])
        
        if not sequences:
            return jsonify({'error': 'No sequences provided'}), 400
        
        if isinstance(sequences, str):
            sequences = [sequences]
            
        print(f"Prediciendo para {len(sequences)} secuencias...")
        
        # Usar la logica existente
        df_result = predictor.predecir(sequences)
        
        if len(df_result) == 0:
             return jsonify([])

        # Convertir dataframe a lista de dicts
        result = df_result.to_dict(orient='records')
        return jsonify(result)

    except Exception as e:
        print(f"Error en prediccion: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/health', methods=['GET'])
def health():
    return jsonify({'status': 'ok', 'service': 'ampclass-api'})

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
