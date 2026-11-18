Python reference
================

The computational API has three entry points. See the :doc:`API guide <api>` for score definitions, input contracts and result-table interpretation.

Evaluation: ITCA
----------------

.. autofunction:: gauge_sc.evaluate

View-specific audit
-------------------

.. autofunction:: gauge_sc.audit

Selection: ITCA+
----------------

.. autofunction:: gauge_sc.select

Plotting
--------

Plot functions return Matplotlib handles. Saving is always an explicit user action.

.. autofunction:: gauge_sc.pl.audit

.. autofunction:: gauge_sc.pl.table

.. autofunction:: gauge_sc.pl.hierarchy

Teaching data
-------------

.. autofunction:: gauge_sc.datasets.audit_example

.. autofunction:: gauge_sc.datasets.load_example

Loading saved results
---------------------

.. autofunction:: gauge_sc.io.load_result

Returned result objects
-----------------------

These classes describe returned objects; users do not need to instantiate or import them from the internal module.

.. autoclass:: gauge_sc._results.EvaluationResult
   :members: save

.. autoclass:: gauge_sc._results.AuditResult
   :members: save, write_to_adata

.. autoclass:: gauge_sc._results.SelectionResult
   :members: save

.. autoclass:: gauge_sc._results.HierarchyEvaluationResult
   :members: save
   :undoc-members:

.. autoclass:: gauge_sc._results.HierarchySelectionResult
   :members: save, write_to_adata
   :undoc-members:
