Installation
============

.. raw:: html

    <style>
        .breadcrumb {
            display: none;
        }
        h1 {
            text-align: center;
            margin-bottom: 15px;
        }
        p.lead.grey-text {
            margin-bottom: 30px;
        }
        .footer-relations {
            border-top: 0px;
        }
        pre {
            background-color: #FFFFFF00;
            border-style: hidden;
            font-size: 87.5% !important;
            margin: 0 !important;
            padding: 0;
        }
        code {
            margin: 0 !important;
            padding: 0 40px 0 40px !important;
        }
    </style>

    <div class="container" id="main-column">
        <div class="text-center">
            <p class="lead grey-text w-responsive mx-auto">
                FlamingPy requires Python 3.8 or newer.
            </p>
            <p class="lead grey-text w-responsive mx-auto mb-6">
                If you currently do not have Python 3 installed, we recommend
                <a href="https://www.anaconda.com/download/">Anaconda for Python 3</a>,
                a distributed version of Python packaged for scientific computation. 
            </p>
            <p class="lead grey-text w-responsive mx-auto mb-6">
                Upon activating a Python environment, run one of the following
                in your choice of CLI: 
            </p>
        </div>

        <ul class="picker nav nav-pills nav-justified mt-6" id="version">
            <li class="nav-item">
                <a class="nav-link active" data-toggle="tab" href="#stable" role="tab">Stable (recommended)</a>
            </li>
            <li class="nav-item">
                <a class="nav-link" data-toggle="tab" href="#dev" role="tab">Source and Testing</a>
            </li>
        </ul>

        <!-- Tab panels -->
        <div class="tab-content pt-0" id="tab-version">
            
            <div class="tab-pane in show active" id="stable" role="tabpanel">
                <pre>
                    <code class="bash">
    # Install the latest released PyPI version of FlamingPy (including
    # dependencies and precompiled C++ binaries), all in a single package.
    python -m pip install flamingpy
                    </code>
                </pre>
            </div>
            <div class="tab-pane slide" id="dev" role="tabpanel">
                <pre>
                    <code class="bash">
    # Download and install the latest source code from GitHub for developers.
    git clone https://github.com/XanaduAI/flamingpy.git
    cd flamingpy
    python -m pip install -r dev-requirements.txt
    # Choose one or some of the following:
    python setup.py develop # Only install Python libraries.
    python setup.py build_cython --inplace # Compile Cython-based backends.
    python setup.py build_cmake --inplace # Compile CMake-based backends.
                    </code>
                </pre>
            </div>
        </div>
    </div>

    <script type="text/javascript">
        $(function(){
            let params = new URLSearchParams(window.location.search);
            let version = params.get("version");

            if (version) {
                $("#version li a").removeClass("active");
                $("#tab-version .tab-pane").removeClass("active");
                $("a[href='#" + version + "']").addClass("active");
                $("#" + version).show();
            };

            $("#version .nav-item a").click(function (e) {
                const old_version = version;
                const new_version = this.hash.substr(1);
                if (old_version != new_version) {
                    $("#" + old_version).hide();
                    $("#" + new_version).show();
                    params.set("version", new_version);
                    const newRelativePathQuery = window.location.pathname + "?" + params.toString();
                    history.pushState(null, "", newRelativePathQuery);
                    version = new_version;
                };
            });

            // Change active navbar element to "Install".
            $(".nav-item.active").removeClass("active");
            $(".nav-item a:contains('Install')").parent().addClass("active");
        });
    </script>
