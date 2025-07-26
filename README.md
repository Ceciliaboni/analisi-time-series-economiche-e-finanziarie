# Analisi di Serie Temporali Economiche e Finanziarie

Questo progetto contiene un'analisi statistica completa di due serie storiche:

1. **Indice dei Prezzi al Consumo (USA, 1985–2024)** – Modello ARIMA
2. **Prezzo Azionario Amazon (2016–2025)** – Modelli GARCH (s-GARCH, GJR-GARCH, T-GARCH, IGARCH)

## Linguaggio e strumenti usati
- 📊 R (con librerie: `forecast`, `rugarch`, `quantmod`, `urca`, `FinTS`, ecc.)
- 📄 Quarto Markdown (`.qmd`) per l'integrazione codice + testo
- 📁 Report esportato in PDF/HTML

## File
- `Programma.qmd` → codice e testo integrati
- `Programma.html` → versione visualizzabile online
- Cartella `R/` → funzioni esterne importate via `source()`

## Output principali
- Stima modelli ARIMA e GARCH
- Diagnostica dei residui
- Previsioni ex-post ed ex-ante
- Confronto modelli tramite test statistici

## Autori
Boni Cecilia & D’Agostino Federica — Università di Firenze, Corso di Laurea in Statistica

📅 Febbraio 2025
