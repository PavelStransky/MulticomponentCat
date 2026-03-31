# Creating multicomponent Schrödinger cat states in a coupled qubit-oscillator system
Accepted as a Letter in the Physical Review A.

### Source codes for generating data and figures

All data for the submitted research paper have been computed using this code.

<ul>
    <li>
        Figure 1: <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/WignerFunction.jl">WignerFunction.jl</a> (Wigner functions and expectation values), <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/ClassicalTrajectory.jl">ClassicalTrajectory.jl</a> (classical trajectories in all subspaces).
    </li>
    <li>
        Figure 2: <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/StrengthFunction.jl">StrengthFunction.jl</a> (strength function and level dynamics).
    </li>
    <li>
        Figure 3: <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/Purity.jl">Purity.jl</a> (Purity), <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/WignerNegativity.jl">WignerNegativity.jl</a>, <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/AnalyseNegativity.py">AnalyseNegativity.jl</a> (the first code calculates a time series of Wigner functions, the second code takes all the elements of the time series and computes the negativity), <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/MeasurementRun.jl">MeasurementRun.jl</a> (purity after the measurement).
    </li>
    <li>
        Figure 4: <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/Periods.nb">Periods.nb</a> (analytical derivations, periods of the classical orbits and classical amplitudes).
    </li>
    <li>
        Figure 5: <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/Measurement.jl">Measurement.jl</a>.
    </li>
    <li>
        Additional code: <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/Sphere.jl">Sphere.jl</a> (the Bloch sphere with a trajectory given by expectation values of J, or other expectation values), <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/CombineFigures.jl">CombineFigures.jl</a> (pastes the spheres into the Wigner function snapshots calculated by <a href="https://github.com/PavelStransky/MulticomponentCat/blob/main/WignerFunction.jl">WignerFunction.jl</a>, resulting in the frames for the animations).
    </li>
</ul>

The remaining source files are either auxiliary files, or were used for additional numerical analyses for investigating in detail the behaviour of the system under other parameter choices.

### Animations
Some animations are published in a <a href="https://www.youtube.com/playlist?list=PLbYyGUIyYU9oRzAqMXd47A2M8ENA4Aumw">YouTube channel</a>.
Animation frames were calculated by this code and put together with <a href="https://www.ffmpeg.org/">FFmpeg</a>.

### Prerequisites
<ul>
    <li>
        <a href="https://qojulia.org/">QuantumOptics.jl</a> Julia package (Rabi model, Wigner functions)
    </li>
    <li>
        <a href="https://docs.sciml.ai/DiffEqDocs">DifferentialEquations.jl</a> Julia package (Classical trajectories)
    </li>
</ul>
    
### Acknowledgements
<ul>
    <li>
        Czech Science Foundation <a href="https://old.starfos.tacr.cz/cs/project/GA22-27431S">25-16056S</a>
    </li>
</ul>

### Final remarks
Feel free to ask the authors for any additional information about the research results or about the code.