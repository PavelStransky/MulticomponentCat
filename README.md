# Creating multicomponent Schrödinger cat states in a coupled qubit-oscillator system
Submitted to Physical Review Letters

### Source codes for generating data and figures

All data for the submitted research paper have been computed using this code.

<ul>
    <li>
        Figure 1: <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/WignerFunction.jl">WignerFunction.jl</a> (Wigner functions), <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/Trajectory.jl">Trajectory.jl</a> (classical trajectory in all subspaces)
    </li>
    <li>
        Figure 2: <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/StrengthFunction.jl">StrengthFunction.jl</a> (strength function and level dynamics)
    </li>
    <li>
        Figure 3: <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/Purity.jl">Purity.jl</a> (Purity), <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/WignerNegativity.jl">WignerNegativity.jl</a>, <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/AnalyseNegativity.py">AnalyseNegativity.jl</a> (the first code calculates a time series of Wigner functions, the second code takes all the time series and computes the negativity)
    </li>
    <li>
        Figure 4: <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/Amplitudes.jl">Amplitudes.jl</a> (amplitudes and Rabi oscillations), <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/Periods.nb">Periods.nb</a> (analytical derivations and periods of classical orbits)
    </li>
</ul>

The remaining source files were used for additional numerical analyses and to investigate in detail the behaviour of the system under other parameter choices.

### Prerequisites and Acknowledgements
<ul>
    <li>
        <a href="https://qojulia.org/">QuantumOptics.jl</a> Julia package
    </li>
    <li>
        Czech Science Foundation <a href="https://old.starfos.tacr.cz/cs/project/GA22-27431S">25-16056S</a>
    </li>
</ul>
    
### Final remarks
Feel free to ask the authors for any additional information about the research results or about the code.