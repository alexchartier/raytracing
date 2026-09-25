# Server privacy

- Anything this project creates or changes on a server must be private to `chartat1`.
- Use `chartat1`'s own account and private home or scratch directories. Do not place project code, jobs, logs, results, or credentials in public or group-accessible trees.
- Set `umask 077`; create directories with mode `0700` and files with mode `0600`. Verify ownership and permissions after staging and after jobs finish.
- Keep scheduler output, temporary files, and copied dependencies under the same private directory. If a server or scheduler cannot enforce this, stop before submitting work.
