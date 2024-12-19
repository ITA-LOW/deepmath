import numpy as np
import matplotlib.pyplot as plt

# --- Parâmetros globais ---
sigma = 10.0
rho = 28.0
beta = 8 / 3
t0, tf = 0, 100  # Tempo inicial e final
h = 0.01       # Passo de tempo
num_passos = int((tf - t0) / h)  # Número total de passos
eta = np.sqrt(beta * (rho - 1))  # Valor de eta
estado_inicial = np.array([rho - 1, eta, eta - 3])  # Condições iniciais


# --- Sistema de Lorenz ---
def sistema_lorenz(estado, sigma, rho, beta):
    x, y, z = estado
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return np.array([dx, dy, dz])


# --- Método de Euler ---
def metodo_euler(estado, h, num_passos):
    x, y, z = estado
    trajetoria = [estado]
    for _ in range(num_passos):
        dx, dy, dz = sistema_lorenz([x, y, z], sigma, rho, beta)
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
        dx1, dy1, dz1 = sistema_lorenz([x, y, z], sigma, rho, beta)
        x_pred = x + h * dx1
        y_pred = y + h * dy1
        z_pred = z + h * dz1

        dx2, dy2, dz2 = sistema_lorenz([x_pred, y_pred, z_pred], sigma, rho, beta)
        x += h / 2 * (dx1 + dx2)
        y += h / 2 * (dy1 + dy2)
        z += h / 2 * (dz1 + dz2)

        trajetoria.append([x, y, z])
    return np.array(trajetoria)


# --- Método RK4 ---
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


# --- Função para calcular resíduos locais ---
def calcular_residuo(trajetoria, h):
    residuos = []
    for i in range(len(trajetoria) - 1):
        derivada_numerica = (trajetoria[i + 1] - trajetoria[i]) / h
        derivada_analitica = sistema_lorenz(trajetoria[i], sigma, rho, beta)
        residuo = np.linalg.norm(derivada_numerica - derivada_analitica)
        residuos.append(residuo)
    return np.array(residuos)


# --- Programa principal ---
trajetoria_euler = metodo_euler(estado_inicial.copy(), h, num_passos)
trajetoria_heun = metodo_heun(estado_inicial.copy(), h, num_passos)
trajetoria_rk4 = runge_kutta_4(estado_inicial.copy(), h, num_passos)

# --- Cálculo dos resíduos ---
residuo_euler = calcular_residuo(trajetoria_euler, h)
residuo_heun = calcular_residuo(trajetoria_heun, h)
residuo_rk4 = calcular_residuo(trajetoria_rk4, h)
tempo = np.linspace(t0, tf, num_passos+1)  # Ajuste para ter 10000 elementos

# --- Visualização dos resíduos ---
plt.figure(figsize=(17, 8))

plt.plot(tempo[:-1], residuo_euler, label='Resíduo - Euler', color='red')
plt.plot(tempo[:-1], residuo_heun, label='Resíduo - Heun', color='blue')
plt.plot(tempo[:-1], residuo_rk4, label='Resíduo - RK4', color='green')

plt.xlabel('Tempo', fontsize=14)
plt.ylabel('Resíduo', fontsize=14)
plt.title('Evolução dos Resíduos dos Métodos Numéricos', fontsize=18)
plt.legend(fontsize=12)
plt.grid()
plt.tight_layout()
plt.show()
