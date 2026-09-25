import unittest
from unittest.mock import patch

from tools.cartman_mcp import policy, remote, server


class CartmanPrivateMcpTests(unittest.TestCase):
    def test_only_chartat1_private_zones_are_exposed(self) -> None:
        self.assertEqual(policy.REMOTE_HOST, "cartman")
        self.assertEqual(set(policy.ZONES), {"repo", "sandbox"})
        self.assertTrue(policy.can_write_path("/homes/chartat1/private_raytracing/runs/test"))
        self.assertFalse(policy.can_write_path("/project/ampere/public/test"))
        self.assertFalse(policy.can_write_path("/homes/chartat1/other/test"))
        names = {definition["name"] for definition in server._tool_definitions()}
        self.assertIn("cartman_qsub_submit", names)
        self.assertNotIn("cartman_bind_mount_ro", names)
        self.assertNotIn("cartman_rsync_from_remote", names)
        self.assertNotIn("cartman_systemctl", names)
        self.assertTrue(remote._bootstrap_remote_shell().startswith("umask 077\n"))

    def test_scheduler_script_stays_owner_only(self) -> None:
        captured = []

        def fake_remote(script: str, timeout_seconds: int = 60) -> remote.RemoteResult:
            captured.append(script)
            return remote.RemoteResult(command=["ssh"], returncode=0, stdout="123.cartman\n", stderr="")

        with patch.object(server, "run_remote_bash", side_effect=fake_remote):
            server._call_tool("cartman_qsub_submit", {
                "zone": "sandbox", "cwd": ".", "script_path": "run/job.sh",
                "script_text": "#!/bin/bash\ntrue", "qsub_args": ["-t", "1-40"],
            })
        self.assertIn("chmod 0700", captured[0])
        self.assertNotIn("0775", captured[0])


if __name__ == "__main__":
    unittest.main()
