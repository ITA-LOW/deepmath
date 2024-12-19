# Atrator de Lorenz - Métodos Numéricos

Este repositório apresenta um código em Python que ilustra a trajetória do **atrator de Lorenz**, resolvendo suas equações diferenciais com três métodos numéricos distintos: Euler, Heun e Runge-Kutta de 4ª ordem (RK4). O atrator de Lorenz é conhecido por seu comportamento dinâmico não-linear e caótico.
## Objetivo

O objetivo deste código é resolver numericamente o sistema de equações diferenciais que descrevem o atrator de Lorenz e visualizar suas trajetórias utilizando três métodos numéricos clássicos. As trajetórias são plotadas em gráficos 3D para demonstrar a diferença entre os métodos utilizados.

## Estrutura do Código

O código é dividido em várias partes principais:

1. **Parâmetros Globais**: São definidos os parâmetros do atrator de Lorenz, como as constantes `sigma`, `rho` e `beta`, além do intervalo de tempo e do passo de tempo utilizado nas simulações.

2. **Métodos Numéricos**: O código implementa três métodos para resolver as equações diferenciais que descrevem o atrator de Lorenz:
   - **Método de Euler**: Método de integração numérica simples, utilizado como base para outros métodos mais sofisticados.
   - **Método de Heun**: Um aprimoramento do método de Euler que utiliza uma correção para obter maior precisão.
   - **Método Runge-Kutta de 4ª ordem (RK4)**: Um método avançado que oferece uma boa precisão para integrações de sistemas dinâmicos.

3. **Visualização**: As trajetórias obtidas pelos três métodos são plotadas em gráficos 3D utilizando a biblioteca `matplotlib`.

## Como Rodar o Código

Para rodar o código, basta ter o Python 3 instalado em sua máquina e as bibliotecas `numpy` e `matplotlib` configuradas.

1. Clone este repositório:

```bash
git clone https://github.com/seu-usuario/atrator-lorenz.git

