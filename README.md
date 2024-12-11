<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GitHub Remote Repo Setup Guide</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 0;
            background-color: #f9f9f9;
            color: #333;
        }
        .container {
            max-width: 800px;
            margin: 0 auto;
            padding: 20px;
            background: #fff;
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
            border-radius: 8px;
        }
        .logo {
            text-align: center;
            margin-bottom: 20px;
        }
        .logo img {
            width: 400px;
            height: auto;
        }
        h1, h2, h3 {
            color: #0056b3;
        }
        pre {
            background: #f4f4f4;
            padding: 10px;
            border-radius: 4px;
            overflow-x: auto;
        }
        code {
            font-family: Consolas, Monaco, 'Courier New', monospace;
        }
        .alert {
            background: #fff4e5;
            border-left: 6px solid #ff9800;
            padding: 10px;
            margin-bottom: 20px;
        }
        ul {
            margin-left: 20px;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="logo">
            <img src="assets/logo.png" alt="Logo">
        </div>

        <h1>A Columnflow-Based Analysis Framework</h1>
        <p><strong>Developed by IPHC and DESY</strong></p>

        <h2>Resources</h2>
        <ul>
            <li><a href="https://github.com/columnflow/columnflow/tree/master" target="_blank">columnflow</a></li>
            <li><a href="https://github.com/riga/law" target="_blank">law</a></li>
            <li><a href="https://github.com/riga/order" target="_blank">order</a></li>
            <li><a href="https://github.com/spotify/luigi" target="_blank">luigi</a></li>
        </ul>

        <h2>Setting Up the Framework (Main Branch)</h2>

        <h3>Step 1: Clone the Repository</h3>
        <p>Clone the repository using the following command:</p>
        <pre><code>git clone --recurse-submodules git@github.com:DesyTau/CPinHToTauTau.git</code></pre>

        <div class="alert">
            <strong>Attention:</strong> If you are setting up the framework on <strong>lxplus</strong> (e.g., on your <code>/afs/cern.ch/user/[first_char]/[cern_user_name]</code>), you might encounter issues with <code>micromamba</code>. 
            <br><br>
            To avoid this, use <code>/eos/user/[first_char]/[cern_user_name]/software</code> as the path for <code>vens</code>, <code>cmssw</code>, and <code>conda</code> folders. The issue is specific to <code>afs</code>, and although the columnflow team is aware of it, a fix is not available (Jacopo, 11/12/24).
        </div>

        <h3>Step 2: Set Up the Framework</h3>
        <p>Once cloned, navigate to the repository directory and set up the framework using:</p>
        <pre><code>cd CPinHToTauTau
source setup.sh [your_setup_name]</code></pre>

        <h2>Setting Up the Framework Used at DESY</h2>
        <p>Follow the same steps outlined above to set up the framework for DESY environments.</p>
    </div>
</body>
</html>



