# Verification-of-Optimal-State-Estimation-for-Optomechanical-Systems

This project aims to verify the conditional state of a quantum optomechanical system through the use of the Wiener filter. The importance of this type of protocol lies in its ability to extract information about the quantum system from the data provided by photodetector measurements. This protocol enables future analyses for real-time manipulation through frequency domain analysis (control theory) or future tests for hypotheses about wave function collapse (fundamental science).

In the code, we use the quantum Wiener filter for both prediction and retrodiction. This combination allows us to access the conditional variance obtained through Bayes' theorem, where the probability of our Gaussian states for position and momentum is given. 

The conditional variance that defines our quantum system is as follows:

![\mathbf{V} = \begin{pmatrix}1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9\end{pmatrix}](https://latex.codecogs.com/png.latex?\mathbf{V}=\begin{pmatrix}V_{11}&V_{12}&;\\;V_{21}&V_{22}\end{pmatrix})

Where ![S+_XX](https://latex.codecogs.com/png.latex?V_{11}(\omega)) ![S+_XX](https://latex.codecogs.com/png.latex?(V_{22}(\omega))) is the conditional variance in position (momentum), and  ![S+_XX](https://latex.codecogs.com/png.latex?V_{12}(\omega)) is the covariance. Based on the information from the power spectral density of the photodetector ![S+_XX](https://latex.codecogs.com/png.latex?S_{XX}(\omega)), we aim to obtain  ![S+_XX](https://latex.codecogs.com/png.latex?V_{11}(\omega)) and  ![S+_XX](https://latex.codecogs.com/png.latex?V_{22}(\omega)). This information can be obtained through the following relationships:

![Equation](https://latex.codecogs.com/svg.image?\int^{\infty}_{-\infty}\frac{1}{2\pi}\big(S_{\overleftarrow{q}\overleftarrow{q}}(\omega)d\omega-S_{\overrightarrow{q}\overrightarrow{q}}(\omega)d\omega\big)\approx&space;2&space;V_{11})

![Equation](https://latex.codecogs.com/svg.image?\int^{\infty}_{-\infty}\frac{1}{2\pi}\big(S_{\overleftarrow{p}\overleftarrow{p}}(\omega)d\omega-S_{\overrightarrow{p}\overrightarrow{p}}(\omega)d\omega\big)\approx&space;2&space;V_{22})

![S+_XX](https://latex.codecogs.com/png.latex?\overrightarrow{q}(\omega)) and ![S+_XX](https://latex.codecogs.com/png.latex?\overleftarrow{q}(\omega)) are the prediction and retrodiction for the position in the frequency domain respectively. This represents our capacity to use the Bayes' theorem in both direction of time to obtain based on the measurement record. To proceed we use the following relation between the Wiener filter and the Power Spectral Density of the outcome in the photodetector : ![S_XXH_q=S_q_q](https://latex.codecogs.com/png.latex?S_{XX}(\omega)\overrightarrow{H}_{q}(\omega)=S_{\overrightarrow{q}\overrightarrow{q}}(\omega)) and![S_XXH_q=S_q_q](https://latex.codecogs.com/png.latex?S_{XX}(\omega)\overleftarrow{H}_{q}(\omega)=S_{\overleftarrow{q}\overleftarrow{q}}(\omega))

Finally, we have an expression where the Wiener filter for prediction and retrodictin is included:

![Equation](https://latex.codecogs.com/svg.image?\int^{\infty}_{-\infty}\frac{1}{2\pi}\big(S_{\overleftarrow{q}\overleftarrow{q}}(\omega)d\omega-S_{\overrightarrow{q}\overrightarrow{q}}(\omega)d\omega\big)S_{XX}\approx&space;2&space;V_{11})

And similarly for the momentum:

![Equation](https://latex.codecogs.com/svg.image?\int^{\infty}_{-\infty}\frac{1}{2\pi}\big(S_{\overleftarrow{p}\overleftarrow{p}}(\omega)d\omega-S_{\overrightarrow{p}\overrightarrow{p}}(\omega)d\omega\big)S_{XX}\approx&space;2&space;V_{22})

Furthermore, the analysis for each power spectral calculation is performed through frequency resolution to evaluate the effect of other frequencies involved in the optomechanical system.
## Plotting Results

The plotting result for the MATLAB code can be seen in the file "Results" from this repository

Please, for a detailed explanation, check this: Santiago-Condori, J. G., Yamamoto, N., Matsumoto, N.(2023). Verification of conditional mechanical squeezing for a mg-scale pendulum near quantum regimes (arXiv:2008.10848).
