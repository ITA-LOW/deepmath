import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

"""
Esse código implementa o método RK4 mais refinado, diminuindo o passo de tempo de 10e-4 para 10e-5. Essa solução foi escolhida para
comparar o quanto as demais soluções implementadas por métodos numéricos diferentes iriam divergir do RK4 mais refinado. 
"""

# --- Parâmetros globais ---
sigma = 10.0
rho = 28.0
beta = 8 / 3
t0, tf = 0, 10  # Tempo inicial e final
h_fino = 0.0001  # Passo de tempo menor para maior precisão
h_grosso = 0.001  # Passo de tempo original
num_passos_fino = int((tf - t0) / h_fino)  # Número total de passos para RK4 refinado
num_passos_grosso = int((tf - t0) / h_grosso)  # Número total de passos para os outros métodos
eta = np.sqrt(beta * (rho - 1))  # Valor de eta
estado_inicial = np.array([rho - 1, eta, eta - 3])  # Condições iniciais

# --- Método de Euler ---
def metodo_euler(estado, h, num_passos):
    x, y, z = estado
    trajetoria = [estado]
    for _ in range(num_passos):
        dx = sigma * (y - x)
        dy = x * (rho - z) - y
        dz = x * y - beta * z
        x += h * dx
        y += h * dy
        z += h * dz
        trajetoria.append([x, y, z])
    return np.array(trajetoria)

# --- Método de Heun ---
def metodo_heun(estado, h, num_passos):
    x, y, z = estado
    trajetoria = [estado]
    for _ in range(num_passos):
        # Predição inicial (Euler)
        dx_euler = sigma * (y - x)
        dy_euler = x * (rho - z) - y
        dz_euler = x * y - beta * z
        x_pred = x + h * dx_euler
        y_pred = y + h * dy_euler
        z_pred = z + h * dz_euler

        # Correção (Heun)
        dx_corr = sigma * (y_pred - x_pred)
        dy_corr = x_pred * (rho - z_pred) - y_pred
        dz_corr = x_pred * y_pred - beta * z_pred
        x += h / 2 * (dx_euler + dx_corr)
        y += h / 2 * (dy_euler + dy_corr)
        z += h / 2 * (dz_euler + dz_corr)
        trajetoria.append([x, y, z])
    return np.array(trajetoria)

# --- Método RK4 ---
def sistema_lorenz(estado, sigma, rho, beta):
    x, y, z = estado
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return np.array([dx, dy, dz])

def runge_kutta_4(estado, h, num_passos):
    trajetoria = [estado]
    for _ in range(num_passos):
        k1 = h * sistema_lorenz(estado, sigma, rho, beta)
        k2 = h * sistema_lorenz(estado + 0.5 * k1, sigma, rho, beta)
        k3 = h * sistema_lorenz(estado + 0.5 * k2, sigma, rho, beta)
        k4 = h * sistema_lorenz(estado + k3, sigma, rho, beta)
        estado = estado + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        trajetoria.append(estado)
    return np.array(trajetoria)

# --- Programa principal ---
# Executar RK4 com maior precisão
trajetoria_rk4_fino = runge_kutta_4(estado_inicial.copy(), h_fino, num_passos_fino)

# Executar os métodos de menor precisão
trajetoria_euler = metodo_euler(estado_inicial.copy(), h_grosso, num_passos_grosso)
trajetoria_heun = metodo_heun(estado_inicial.copy(), h_grosso, num_passos_grosso)
trajetoria_rk4_grosso = runge_kutta_4(estado_inicial.copy(), h_grosso, num_passos_grosso)

# --- Interpolação das trajetórias ---
# Tempo para RK4 refinado
t_fino = np.linspace(t0, tf, num_passos_fino + 1)

# Interpolação das trajetórias dos métodos de menor precisão
interp_euler = interp1d(np.linspace(t0, tf, num_passos_grosso + 1), trajetoria_euler, axis=0)
interp_heun = interp1d(np.linspace(t0, tf, num_passos_grosso + 1), trajetoria_heun, axis=0)
interp_rk4 = interp1d(np.linspace(t0, tf, num_passos_grosso + 1), trajetoria_rk4_grosso, axis=0)

# Trajetórias interpoladas
trajetoria_euler_interp = interp_euler(t_fino)
trajetoria_heun_interp = interp_heun(t_fino)
trajetoria_rk4_interp = interp_rk4(t_fino)

# --- Calcular erro global para cada método ---
def calcular_erro_global(trajetoria_ref, trajetoria):
    # Usando a norma L2 para calcular o erro global
    erro = np.linalg.norm(trajetoria_ref - trajetoria, axis=1)
    return np.mean(erro)

# Erros globais
erro_euler = calcular_erro_global(trajetoria_rk4_fino, trajetoria_euler_interp)
erro_heun = calcular_erro_global(trajetoria_rk4_fino, trajetoria_heun_interp)
erro_rk4 = calcular_erro_global(trajetoria_rk4_fino, trajetoria_rk4_interp)

# Exibir os erros
print(f"Erro Euler: {erro_euler:.6f}")
print(f"Erro Heun: {erro_heun:.6f}")
print(f"Erro RK4 grosseiro: {erro_rk4:.6f}")
